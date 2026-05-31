/**********************************************************************\
| A Quick tree Drawer - Reads a tree in parenthetical notation (TNT    |
| format, possibly with branch labels), and draws the tree diagram     |
|                by Pablo A. Goloboff (2026)                           |
|         =============================================                |
|               Intended for pedagogic use                             |
| Feel free to use or modify; if you do, please acknowledge authorship |
|         =============================================                |
|       The tree-printing algorithm is taken almost verbatim           |
|           from Goloboff (2022 ISBN 9780367420277)                    |
|   https://www.lillo.org.ar/phylogeny/eduscripts/treeplotting.pic     |
\**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h> 
#include <fcntl.h>
#include <io.h>

#define MAXT 5000
FILE * tfile , * namfile ; 
int curnodin , atnt , nt , anc[2*MAXT] , lef[2*MAXT] , rig[2*MAXT] , sis[2*MAXT] , trylist[2*MAXT] , atlin[2*MAXT] , linsdone ; 
#ifdef POORMAN
int rho = ',' , el = '`' , vert = '|' , hor = '-' , side = '|' , cross = '+' , up = '-' , hook = '|' ; 
#else 
int rho = 218 , el = 192 , vert = 179 , hor = 196 , side = 195 , cross = 197 , up = 193 , hook = 180 ;
#endif
char names[ MAXT ] [ 150 ] , msg[ 2*MAXT ] [ 50 ] , tmpstr[ 150 ];
int * retrop , * retrolim , retrolist[ MAXT ] , leglen[2*MAXT] , cladesize[2*MAXT] , usernames = 0 ; 

/*** SOME SUPPORT FUNCTIONS ******************/

void errout ( char * txt ) {
    fprintf ( stderr , "Can't run:\n    %s\n" , txt ) ;
    getc ( stdin ) ; 
    exit ( 0 ) ; }

int listabove( int * list , int from ) {
    int * done , * todo , at , n ;
    * ( done = list ) = from ;
    todo = done + 1 ;
    while ( done < todo ) {
       at = * done ++ ;
       if ( at >= nt ) 
          for ( n = lef[ at ] ; n >= 0 ; n = sis [ n ] )  
              if ( n >= nt ) * todo ++ = n ; }
    return done - list ; }

/***  READING THE TREE FROM FILE ************/

void add ( int what , int where ) {
    anc [ what ] = where ; 
    if ( lef [ where ] < 0 ) lef[ where ] = rig[ where ] = what ;
    else {
        sis [ rig [ where ] ] = what ;
        rig[ where ] = what ; }}

int stringto( char * str , int matchnames ) {
   char * cp , c ;
   int i ; 
   if ( matchnames ) cp = tmpstr - 1 ; 
   else cp = str - 1 ;  
   c = 'a' ; 
   while ( !isspace ( c ) && c != '=' && c != '(' && c != ')' && c != '*' && c != ';' && !feof ( tfile ) ) {
       * ++ cp = c = getc ( tfile ) ;
       if ( matchnames && !isspace ( c ) && c != '=' && c != '(' && c != ')' && c != '*' && c != ';' && !isdigit( c ) ) errout ( "When giving a list of names,\n    tree must consist of numbers only" ) ; }
   * cp = '\0' ;
   if ( matchnames ) {
       i = atoi ( tmpstr ) ;
       if ( i >= nt ) errout ( "No taxon %i in name file" ) ; }
   if ( feof ( tfile ) ) errout ( "Unexpected end-of-file" ) ; 
   if ( !isspace ( c ) ) ungetc( c , tfile ) ;
   if ( matchnames ) return i ;
   return 0 ; }

void treein ( int at ) {
    char c ;
    int i , n ;
    while ( 1 ) {
        fscanf ( tfile , " %c" , &c ) ;
        if ( c == ')' ) {
            fscanf ( tfile , " %c" , &c ) ;
            if ( c == '=' ) stringto( msg[ at ] , 0 ) ; 
            else ungetc ( c , tfile ) ; 
            return ; }
        if ( c == '(' ) {
            if ( ++ curnodin >= 2*nt-1 ) errout ( "Too many parentheses" ) ;
            add ( curnodin , at ) ;
            treein ( curnodin ) ; }
        else {
            if ( c == '=' ) { i = atnt - 1 ; fscanf ( tfile , "%s" , &msg[ i ] ) ; }
            else {
               ungetc ( c , tfile ) ;
               i = atnt ;
               n = stringto ( names[ atnt ++ ] , usernames ) ;
               if ( usernames ) i = n ; 
               add ( i , at ) ; }}}}  
    
void open_tree_file ( char * fnam ) {
    int i ;
    char command[150] ; 
    if ( ( tfile = fopen ( fnam , "rb" ) ) == NULL ) errout ( "Tree file" ) ;
    fscanf ( tfile , " %s" , &command ) ;
    if ( strcmp ( command , "tread" ) ) errout ( "File does not start with \"tread\"" ) ;
    i = ' ' ; while ( isspace ( i ) ) i = getc ( tfile ) ;
    if ( i == 39 ) {
        i = ' ' ; while ( i != 39 ) i = getc ( tfile ) ;
        i = ' ' ; while ( isspace ( i ) ) i = getc ( tfile ) ; }
    if ( i == '(' ) ungetc ( i , tfile ) ;
    else errout ( "Expected parenthesis right after \"tread\"" ) ; }
    
char read_tree ( void ) {
    int i ;
    i = ' ' ; while ( isspace ( i ) ) i = getc ( tfile ) ; 
    if ( i != '(' ) errout ( "Parentheses" ) ;
    atnt = 0 ;
    for ( i = 0 ; i < 2*MAXT ; ++ i ) msg[ i ] [ 0 ] = 0 ; // make sure all legends blank by default!
    for ( i = 0 ; i < 2*nt-1 ; ++ i ) anc[ i ] = lef[ i ] = rig[ i ] = sis[ i ] = -1 ;
    treein ( curnodin = nt ) ;
    anc[ nt ] = ++ curnodin ;
    i = ' ' ; while ( isspace ( i ) ) i = getc ( tfile ) ; 
    return i ; }

void gettaxnum ( char * fnam ) {
    int i ;
    char command[150] ; 
    if ( ( tfile = fopen ( fnam , "rb" ) ) == NULL ) errout ( "Tree file" ) ;
    fscanf ( tfile , " %s" , &command ) ;
    if ( strcmp ( command , "tread" ) ) errout ( "File does not start with \"tread\"" ) ;
    i = ' ' ; while ( isspace ( i ) ) i = getc ( tfile ) ;
    if ( i == 39 ) {
        i = ' ' ; while ( i != 39 ) i = getc ( tfile ) ;
        i = ' ' ; while ( isspace ( i ) ) i = getc ( tfile ) ; }
    if ( i != '(' ) errout ( "Expected parenthesis right after \"tread\"" ) ; 
    nt = 0 ;      
    while ( 1 ) {
        i = getc( tfile ) ;
        if ( isspace ( i ) ) continue ; 
        if ( i == ';' || i == '*' ) break ; 
        if ( i == '(' || i == ')' ) continue ;
        if ( i == '=' ) { stringto ( command , 0 ) ; continue ; }
        ungetc ( i , tfile ) ;
        stringto( command , 0 ) ; 
        ++ nt ; }
    fclose ( tfile ) ; }

/***  PREPARING THE TREE FOR PLOTTING **********/

void initlines ( int which )   // Find out at which line every terminal is printed 
{ int i ; 
  if ( which >= nt ) {
     for ( i = lef[ which ] ; i >= 0 ; i = sis [ i ] ) initlines ( i ) ;  
     return ; }
  atlin [ which ] = ++ linsdone ; }     

void calculines ( int which )   // Find out at which line every internal node is printed 
{ int i ; 
  for ( i = lef[ which ] ; i >= 0 ; i = sis [ i ] ) 
      if ( i >= nt ) calculines ( i ) ;  
  atlin [ which ] = ( 1 + atlin [ lef [ which ] ] +  atlin [ rig [ which ] ] ) / 2 ; }     

void prepare_tree ( void ) 
{ int i , j , k , n , q , switched ; 
  for ( i = 0 ; i < nt ; ++ i ) cladesize[ i ] = 1 ; 
  n = listabove( trylist , nt ) ; 
  for ( k = n ; k -- ; ) {
    i = trylist[ k ] ; 
    for ( cladesize[ i ] = 0 , j = lef[ i ] ; j >= 0 ; j = sis [ j ] ) 
         cladesize [ i ] += cladesize[ j ] ; }
  for ( i = nt ; i < anc[ nt ] ; ++ i ) {
      for ( k = 0 , n = lef[ i ] ; n >= 0 ; n = sis [ n ] , k ++ ) 
          trylist[ k ] = n ;
      for ( n = switched = 0 ; n < k - 1 ; ++ n ) 
        for ( j = n + 1 ; j < k ; ++ j ) {  // a bubble sort, but this is expected to be small numbers... should improve for large trees
          if ( cladesize[ trylist[ j ] ] < cladesize[ trylist[ n ] ] ) {
             switched = 1 ; 
             q = trylist[ n ] ;
             trylist[ n ] = trylist[ j ] ;
             trylist[ j ] = q ; }}
      if ( switched ) {
          lef[ i ] = trylist[ 0 ] ;
          rig [ i ] = trylist[ k - 1 ] ;
          sis[ trylist[ k - 1 ] ] = -1 ; 
          for ( n = 0 ; n < k - 1 ; n ++ )
             sis [ trylist[ n ] ] = trylist[ n + 1 ] ; }}
  for ( i = 0 ; i < anc[ nt ] ; ++ i ) 
      leglen[ i ] = strlen ( msg[ i ] ) ; 
  retrolim = retrolist ; 
  retrolim += nt ;
  linsdone = 0 ; 
  initlines ( nt ) ; 
  calculines ( nt ) ; 
  linsdone = 0 ; }

/**  PLOTTING THE TREE ************/

void setretrop ( int which ) 
{ int i ; 
  retrop = retrolim ; 
  i = which ; 
  while ( i != nt ) {
      i = anc[ i ] ; 
      * --  retrop = i ; }}

void plot ( int which ) 
{ int i , j ; 
  int joint ; 
  if ( which >= nt ) {
      for ( i = lef[ which ] ; i >= 0 ; i = sis [ i ] ) plot ( i ) ;  
      return ; }
  setretrop ( which ) ; 
  linsdone ++ ; 
  while ( 1 ) {
      j = * retrop ++ ;
      if ( linsdone == atlin [ j ] ) {
             printf( "%c%c" , hor , hor ) ;
             if ( msg[ j ] [ 0 ] ) printf ( "%s" , msg [ j ] ) ; }
      else { printf ( "  " ) ; 
             if ( msg[ j ] [ 0 ] ) for ( i = 0 ; i < leglen [ j ] ; ++ i ) printf( " " ) ; }
      joint = ' ' ; 
      if ( linsdone <= atlin [ rig [ j ] ] ) {
          if ( linsdone >= atlin [ lef [ j ] ] ) joint = vert ; 
          if ( linsdone == atlin [ j ] ) joint = hook ; 
          for ( i = lef [ j ] ; i >= 0 ; i = sis [ i ] ) 
              if ( atlin [ i ] == linsdone ) {
                  if ( i == lef [ j ] ) joint = rho ; 
                  else if ( i == rig [ j ] ) 
                          if ( cladesize[ j ] == 2 ) joint = up ; 
                          else joint = el ; 
                  else if ( linsdone == atlin [ j ] ) joint = cross ; 
                       else joint = side ; }}
      printf ( "%c" , joint ) ;       
      if ( retrop == retrolim ) break ; }
  printf( "%c" , hor ) ;     
  if ( msg[ which ] [ 0 ] ) printf ( "%c%s" , hor , msg [ which ] ) ; 
  else printf ( "%c" , hor ) ; 
  printf( " %s\n" , names[ which ] ) ; }  

void gettaxnames( char * fnam ) {
    int i ; 
    if ( ( namfile = fopen ( fnam , "rb" ) ) == NULL ) errout ( "Can't open name file" ) ;
    usernames = 1 ; 
    for ( i = 0 ; i < nt && !feof( namfile ) ; ++ i ) 
        fscanf ( namfile , " %s" , &names[ i ] ) ;
    if ( i < nt ) errout ( "Not enough taxon names in file" ) ; }

void main ( int argc , char ** argv ) {
    char endis = '*' ;
    int treen = 0 ; 
    if ( argc == 1 ) errout ( "Give file name (parenthetical tree, TNT format)\n    If tree consists of only numbers, can give\n    list of taxon names in a separate file (2nd argument)" ) ;
    gettaxnum ( argv[ 1 ] ) ;
    open_tree_file ( argv[ 1 ] ) ;
    if ( argc > 2 ) gettaxnames ( argv[ 2 ] ) ; 
    while ( endis == '*' ) {
        endis = read_tree () ; 
        prepare_tree () ; 
        plot ( nt ) ;
        if ( isatty ( fileno ( stdout ) ) ) {
           if ( endis == '*' ) printf ( "\nPress <enter> to show tree %i\n" , ++ treen ) ;  
           else printf ( "\nPress <enter> to exit\n" ) ;  
           getc ( stdin ) ; }}
    fclose ( tfile ) ; }

