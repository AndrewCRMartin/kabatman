BOOL BuildWhere(char *buffer)
;
BOOL CheckForSetOper(char *word, BOOL *error)
;
char *HandleWhereSubClause(char *buffer, char *word, BOOL *error,
                           int maxlength)
;
BOOL SetComparison(WHERE *p, char *word)
;
void SetWhereData(WHERE *wh, char *word)
;
void ClearWhere(void)
;
