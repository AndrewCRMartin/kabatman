BOOL ExecuteSearch(char *filename)
;
BOOL HandleLogical(WHERE *wh, int *StackDepth)
;
BOOL HandleMatch(WHERE *wh, int *StackDepth)
;
BOOL DoStrTest(char *text, int comparison, char *subtext, BOOL fuzzy)
;
BOOL fuzzystrstr(char *text, char *subtext)
;
BOOL DoIntTest(int number, int comparison, int testnum)
;
BOOL DoCharTest(char ch, int comparison, char testch)
;
BOOL DoBoolTest(BOOL condition, int comparison, char *test)
;
BOOL IsComplete(DATA *d)
;
void DisplaySearch(FILE *fp, int StackDepth)
;
void FillLoop(char *loopname, DATA *d, char *loop)
;
char GetResidue(DATA *d, char *resid)
;
BOOL FindCanonical(DATA *d, char *LoopID, char *class)
;
void RemoveDupes(int StackDepth)
;
BOOL TooSimilar(DATA *d, DATA *e, REAL Cutoff)
;
void GetSubgroup(DATA *d, char *chain, char *subgroup)
;
void DoGetSubgroup(char *sequence, char *class, char *subgroup)
;
void WriteAsPIR(FILE *fp, DATA *d)
;
