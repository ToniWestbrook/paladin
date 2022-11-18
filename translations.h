#ifndef TRANSLATIONS_H_
#define TRANSLATIONS_H_

#define TRANSLATION_COUNT 25

extern unsigned char codon_aa_hash[][64];
extern const char * translationNames[];

int * convertTransArgs(const char * passArguments);
void renderTranslations();

#endif /* TRANSLATIONS_H_ */
