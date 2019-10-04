int main(int argc, char **argv)
;
void DisplayCopyright(BOOL DoHash)
;
BOOL ParseCmdLine(int argc, char **argv, BOOL *ForceRead)
;
BOOL ReadStoredData(char *filename)
;
char **ReadSpecialNumbering(char *buffer)
;
BOOL StoreKabatData(char *filename)
;
BOOL ReadKabatData(char *FoF)
;
BOOL StoreKabatInData(DATA **pData, KABATENTRY Kabat, char chain,
                      char *source, BOOL GotInsert)
;
BOOL StoreDataInData(DATA **pData, DATA *p, char chain, char *source)
;
void AddKabatToData(DATA *Data, KABATENTRY Kabat, char chain, 
                    BOOL DoInsert)
;
void AddDataToData(DATA *Data, DATA *Extra, char chain)
;
DATA *StoreHAndMatchL(DATA *Data, KABATENTRY KabatH, DATA *KabL[], 
                      int NLFile, char *source, BOOL GotInsert)
;
DATA *StoreUnmatchedL(DATA *Data, DATA *KabL[], int NLFile, char *source)
;
void CopyKabatToData(DATA *p, KABATENTRY Kabat, char chain, 
                     BOOL GotInsert)
;
void CopyDataToData(DATA *p, DATA *extra, char chain)
;
BOOL ReadLFiles(FILE *fpL[], int NLFile, DATA *KabatLData[], char *source,
                int *LCClass)
;
void GetSource(char *filename, char *source)
;
void CommandLoop(void)
;
BOOL ReadChothiaData(char *filename)
;
BOOL RefCheck(char *inref1, char *inref2)
;
BOOL NameCheck(char *instr1, char *instr2)
;
char *GetFullStopWord(char *buffer, char *word)
;
char *EatInitials(char *buffer)
;
void BuildFrom(char *buffer)
;
void HandleSetCommand(char *buffer)
;
