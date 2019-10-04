int ReadNextKabatEntry(FILE *fp, KABATENTRY *Kabat, BOOL *insert,
                       BOOL OldFormat)
;
void ClearKabatEntry(KABATENTRY *KabatEntry)
;
int GetKabatOffset(char **table, char *label)
;
char **BuildKabatNumbering(KABATENTRY Kabat, BOOL OldFormat)
;
void CheckKabatNumbering(char **KabatIndex)
;
char *CreateKabatNumber(char *resnum, int offset)
;
int ReadKabatEntry(FILE *fp, char *buffer, int bufflen, 
                   KABATENTRY *KabatEntry, BOOL *insert)
;
int ReadOldKabatEntry(FILE *fp, char *buffer, int bufflen, 
                      KABATENTRY *KabatEntry, BOOL *insert)
;
