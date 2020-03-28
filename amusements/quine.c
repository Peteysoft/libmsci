#include <stdio.h>
int main() {char listing[200]="#include <stdio.h>%cint main() {char listing[200]=%c%s%c; char delim=10; char delim2=34; printf(listing, delim, delim2, listing, delim2);}"; char delim=10; char delim2=34; printf(listing, delim, delim2, listing, delim2);}
