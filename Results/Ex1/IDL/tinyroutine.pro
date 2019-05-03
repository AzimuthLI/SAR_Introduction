PRO tinyRoutine 
; Create a string variable 
myString = 'This is a tiny IDL routine' 
PRINT, mystring 
; Create some other variables 
myNumber = 4 
myResult = STRING(myNumber * !PI) 
; Display the myResult variable 
void = DIALOG_MESSAGE('Result: '+myResult) 
END 