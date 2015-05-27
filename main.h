#ifndef MAIN_H_
#define MAIN_H_

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1.0"
#endif

// Render usage and version details
int renderUsage();
int renderVersion();

// CLEAN
int bwa_fa2pac(int argc, char *argv[]);
int main_shm(int argc, char *argv[]);
int main_pemerge(int argc, char *argv[]);

#endif /* MAIN_H_ */
