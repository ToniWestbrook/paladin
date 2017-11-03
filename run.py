#!/usr/bin/python
"""Wrapper script for running PALADIN and returning the results."""

import os
import json
import uuid
import boto3
import shutil
import logging
import argparse
import subprocess


def run_cmds(commands, retry=0, catchExcept=False):
    """Run commands and write out the log, combining STDOUT & STDERR."""
    logging.info("Commands:")
    logging.info(' '.join(commands))
    p = subprocess.Popen(commands,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    stdout, stderr = p.communicate()
    exitcode = p.wait()
    if stdout:
        logging.info("Standard output of subprocess:")
        for line in stdout.split('\n'):
            logging.info(line)
    if stderr:
        logging.info("Standard error of subprocess:")
        for line in stderr.split('\n'):
            logging.info(line)

    # Check the exit code
    if exitcode != 0 and retry > 0:
        msg = "Exit code {}, retrying {} more times".format(exitcode, retry)
        logging.info(msg)
        run_cmds(commands, retry=retry - 1)
    elif exitcode != 0 and catchExcept:
        msg = "Exit code was {}, but we will continue anyway"
        logging.info(msg.format(exitcode))
    else:
        assert exitcode == 0, "Exit code {}".format(exitcode)


def run_paladin(input_str,
                db_fp,
                db_url,
                output_folder,
                threads=16,
                temp_folder='/share/temp'):
    """Align a set of reads against a reference database with Paladin."""

    # Use the read prefix to name the output and temporary files
    read_prefix = input_str.split('/')[-1]

    # Check to see if the output already exists, if so, skip this sample
    output_fp = output_folder.rstrip('/') + '/' + read_prefix + '.json.gz'
    if output_fp.startswith('s3://'):
        # Check S3
        logging.info("Making sure that the output path doesn't already exist")
        bucket = output_fp[5:].split('/')[0]
        prefix = '/'.join(output_fp[5:].split('/')[1:])
        client = boto3.client('s3')
        results = client.list_objects(Bucket=bucket, Prefix=prefix)
        if 'Contents' in results:
            msg = "Output already exists, skipping ({})"
            logging.info(msg.format(output_fp))
            return
    else:
        # Check local filesystem
        if os.path.exists(output_fp):
            msg = "Output already exists, skipping ({})"
            logging.info(msg.format(output_fp))
            return

    # Get the reads
    read_fp = get_reads_from_url(input_str, temp_folder)

    # Align the reads against the reference database
    logging.info("Aligning reads")
    output_prefix = os.path.join(temp_folder, read_prefix)
    run_cmds(["/usr/bin/paladin/paladin",
              "align",
              "-t", str(threads),          # Threads
              "-o", output_prefix,         # Output path
              "-u", "0",                   # Don't contact uniprot.org
              db_fp,                       # Database prefix
              read_fp])                  # FASTQ input
    # Output path
    output_fp = os.path.join(temp_folder, read_prefix + "_uniprot.tsv")
    assert os.path.exists(output_fp)

    # Parse the alignment to get the abundance summary statistics
    logging.info("Parsing the output ({})".format(output_fp))
    paladin_results = parse_tsv(output_fp)

    # Clean up the output and FASTQ
    os.unlink(output_fp)
    # Don't delete local files
    if any([input_str.startswith("sra://"),
            input_str.startswith("s3://"),
            input_str.startswith("ftp://")]):
        os.unlink(read_fp)

    # Read in the logs
    logging.info("Reading in the logs")
    logs = open(log_fp, 'rt').readlines()

    # Make an object with all of the results
    out = {
        "input_path": input_str,
        "input": read_prefix,
        "output_folder": output_folder,
        "logs": logs,
        "ref_db": db_fp,
        "results": paladin_results
    }

    # Write out the final results JSON and write them to the output folder
    return_results(out, read_prefix, output_folder, temp_folder)


def parse_tsv(fp):
    """Parse the Paladin output (in TSV format)."""
    dat = []
    with open(fp, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            line = line.rstrip("\n").split("\t")
            # Skip empty lines
            if len(line) == 1:
                continue
            assert len(line) == len(header)
            dat.append(dict(zip(header, line)))
    logging.info("Read in {} lines from {}".format(len(dat), fp))
    return dat


def get_sra(accession, temp_folder):
    """Get the FASTQ for an SRA accession via ENA."""
    local_path = os.path.join(temp_folder, accession + ".fastq")
    # Download from ENA via FTP
    # See https://www.ebi.ac.uk/ena/browse/read-download for URL format
    url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq"
    folder1 = accession[:6]
    url = "{}/{}".format(url, folder1)
    if len(accession) > 9:
        if len(accession) == 10:
            folder2 = "00" + accession[-1]
        elif len(accession) == 11:
            folder2 = "0" + accession[-2:]
        elif len(accession) == 12:
            folder2 = accession[-3:]
        else:
            logging.info("This accession is too long: " + accession)
            assert len(accession) <= 12
        url = "{}/{}".format(url, folder2)
    # Add the accession to the URL
    url = "{}/{}/{}".format(url, accession, accession)
    logging.info("Base info for downloading from ENA: " + url)
    # There are three possible file endings
    file_endings = ["_1.fastq.gz", "_2.fastq.gz", ".fastq.gz"]
    # Try to download each file
    for end in file_endings:
        run_cmds(["curl",
                  "-o", os.path.join(temp_folder, accession + end),
                  url + end], catchExcept=True)
    # If none of those URLs downloaded, fall back to trying NCBI
    if any([os.path.exists("{}/{}{}".format(temp_folder, accession, end))
            for end in file_endings]):
        # Combine them all into a single file
        logging.info("Combining into a single FASTQ file")
        with open(local_path, "wt") as fo:
            cmd = "gunzip -c {}/{}*fastq.gz".format(temp_folder, accession)
            gunzip = subprocess.Popen(cmd, shell=True, stdout=fo)
            gunzip.wait()

        # Clean up the temporary files
        logging.info("Cleaning up temporary files")
        for end in file_endings:
            fp = "{}/{}{}".format(temp_folder, accession, end)
            if os.path.exists(fp):
                os.unlink(fp)
    else:
        logging.info("No files found on ENA, trying SRA")
        run_cmds(["fastq-dump", "--outdir", temp_folder, accession])

        # Check to see if the file was downloaded
        msg = "File could not be downloaded from SRA: {}".format(accession)
        assert os.path.exists(local_path), msg

    # Return the path to the file
    logging.info("Done fetching " + accession)
    return local_path


def get_reads_from_url(input_str, temp_folder):
    """Get a set of reads from a URL -- return the downloaded filepath and file prefix."""
    logging.info("Getting reads from {}".format(input_str))

    if input_str.startswith(('s3://', 'sra://', 'ftp://')) is False:
        logging.info("Path does not start with s3://, sra://, or ftp://")

        # Check that the path exists locally
        assert os.path.exists(input_str), "Path does not exist locally"
        logging.info("{} is a valid local path".format(input_str))
        # Return the input string as the valid local path
        return input_str

    filename = input_str.split('/')[-1]
    local_path = os.path.join(temp_folder, filename)

    logging.info("Filename: " + filename)
    logging.info("Local path: " + local_path)

    # Get files from AWS S3
    if input_str.startswith('s3://'):
        logging.info("Getting reads from S3")
        run_cmds(['aws', 's3', 'cp', '--quiet', '--sse', 'AES256', input_str, temp_folder])
        return local_path

    # Get files from an FTP server
    elif input_str.startswith('ftp://'):
        logging.info("Getting reads from FTP")
        run_cmds(['wget', '-P', temp_folder, input_str])
        return local_path

    # Get files from SRA
    elif input_str.startswith('sra://'):
        accession = filename
        logging.info("Getting reads from SRA: " + accession)
        local_path = get_sra(accession, temp_folder)

        return local_path

    else:
        raise Exception("Did not recognize prefix for reads: " + input_str)


def get_reference_database(ref_db, temp_folder):
    """Get a reference database."""

    # Get files from AWS S3
    if ref_db.startswith('s3://'):
        logging.info("Getting reference database from S3: " + ref_db)

        # Save the database to a local path with a random string prefix
        # This avoids collision between multiple running processes
        rand_string = uuid.uuid4()
        local_folder = os.path.join(temp_folder, "{}.db/".format(rand_string))
        os.mkdir(local_folder)

        logging.info("Saving database to " + local_folder)
        run_cmds(['aws', 's3', 'sync', '--quiet', '--sse', 'AES256',
                  ref_db, local_folder])

        # If the database was downloaded from S3, delete it when finished
        delete_db_when_finished = True

        # Get the prefix for the database
        for fp in os.listdir(local_folder):
            if fp.endswith(".pro"):
                prefix = fp[:-4]
                local_fp = os.path.join(local_folder, prefix)
        logging.info("Prefix for reference database is {}".format(local_fp))

        return local_fp, delete_db_when_finished

    else:
        # Treat the input as a local path
        logging.info("Getting reference database from local path: " + ref_db)
        assert os.path.exists(ref_db)

        # Get the prefix for the database
        local_fp = None
        for fp in os.listdir(ref_db):
            if fp.endswith(".pro"):
                prefix = fp[:-4]
                local_fp = os.path.join(ref_db, prefix)
        msg = "No Paladin database could be found in " + ref_db
        assert local_fp is not None, msg
        logging.info("Prefix for reference database is {}".format(ref_db))

        # Don't delete this database when finished
        delete_db_when_finished = False

        return local_fp, delete_db_when_finished


def return_results(out, read_prefix, output_folder, temp_folder):
    """Write out the results as JSON and copy to the output folder."""
    # Make a temporary file
    temp_fp = os.path.join(temp_folder, read_prefix + '.json')
    with open(temp_fp, 'wt') as fo:
        json.dump(out, fo)
    # Compress the output
    run_cmds(['gzip', temp_fp])
    temp_fp = temp_fp + '.gz'

    if output_folder.startswith('s3://'):
        # Copy to S3
        run_cmds(['aws', 's3', 'cp', '--quiet', '--sse', 'AES256',
                  temp_fp, output_folder])
    else:
        # Copy to local folder
        run_cmds(['mv', temp_fp, output_folder])


def make_scratch_space(scratch_size, temp_folder):
    """Create scratch space using a ramdisk."""
    run_cmds(['mount', '-t', 'tmpfs', '-o', 'size={}g'.format(scratch_size),
              'tmpfs', temp_folder])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Align a set of reads against a reference database with DIAMOND, and save the results.
    """)

    parser.add_argument("--input",
                        type=str,
                        help="""Location for input file(s). Comma-separated.
                                (Supported: sra://, s3://, or ftp://).""")
    parser.add_argument("--ref-db",
                        type=str,
                        help="""Folder containing reference database.
                                (Supported: s3://, ftp://, or local path).""")
    parser.add_argument("--output-folder",
                        type=str,
                        help="""Folder to place results.
                                (Supported: s3://, or local path).""")
    parser.add_argument("--scratch-size",
                        type=int,
                        default=None,
                        help="If specified, create a ramdisk of this size (Gb).")
    parser.add_argument("--threads",
                        type=int,
                        default=16,
                        help="Number of threads to use aligning.")
    parser.add_argument("--temp-folder",
                        type=str,
                        default='/share',
                        help="Folder used for temporary files (and ramdisk, if specified).")

    args = parser.parse_args()

    # Set up logging
    log_fp = 'log.txt'
    logFormatter = logging.Formatter('%(asctime)s %(levelname)-8s [run.py] %(message)s')
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    # Write to file
    fileHandler = logging.FileHandler(log_fp)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    # Also write to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    # Set up the scratch space
    if args.scratch_size is not None:
        logging.info("Setting up scratch space ({}Gb)".format(args.scratch_size))
        make_scratch_space(args.scratch_size, args.temp_folder)

    # Get the reference database
    db_fp, delete_db_when_finished = get_reference_database(args.ref_db,
                                                            args.temp_folder)
    logging.info("Reference database: " + db_fp)

    # Align each of the inputs and calculate the overall abundance
    for input_str in args.input.split(','):
        logging.info("Processing input argument: " + input_str)
        run_paladin(input_str,              # ID for single sample to process
                    db_fp,                  # Local path to DB
                    args.ref_db,            # URL of ref DB, used for logging
                    args.output_folder,     # Place to put results
                    threads=args.threads,
                    temp_folder=args.temp_folder)

    # Delete the reference database
    if delete_db_when_finished:
        logging.info("Deleting reference database: {}".format(db_fp))
        shutil.rmtree(db_fp)

    # Stop logging
    logging.info("Done")
    logging.shutdown()
