#!/bin/bash

here=`pwd`

echo "
<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 3.2 Final//EN\">
<html>
 <head>
  <title>Index of ${here}</title>
 </head>
 <body>
<h1>Index of ${here}</h1>
  <table>
   <tr><th valign=\"top\"><img src=\"/icons/blank.gif\" alt=\"[ICO]\"></th><th><a href=\"?C=N;O=D\">Name</a></th><th><a href=\"?C=M;O=A\">Last modified</a></th><th><a href=\"?C=S;O=A\">Size</a></th><th><a href=\"?C=D;O=A\">Description</a></th></tr>
   <tr><th colspan=\"5\"><hr></th></tr>
<tr><td valign=\"top\"><img src=\"/icons/back.gif\" alt=\"[PARENTDIR]\"></td><td><a href=\"..\">Parent Directory</a></td><td>&nbsp;</td><td align=\"right\"> - </td><td>&nbsp;</td></tr>"

rm listing$$.tmp 2>/dev/null
find . -maxdepth 1 -type f -printf "%A@ %f| %AY-%Am-%Ad %.5AT| %s\n" > listing$$.tmp
find . -maxdepth 1 -type d -printf "%A@ %f| %AY-%Am-%Ad %.5AT| - \n" >> listing$$.tmp

sort -n -k 1 listing$$.tmp | cut -d ' ' -f2- | awk -F'|' '{printf("<tr><td valign=\"top\"><img src=\"/icons/folder.gif\" alt=\"[DIR]\"></td><td><a href=\"%s\">%s</a></td><td align=\"right\"> %s </td><td align=\"right\"> %s </td><td>&nbsp;</td></tr>\n",$1,$1,$2,$3)}'

rm listing$$.tmp

echo " <tr><th colspan=\"5\"><hr></th></tr>
</table>
<address>Apache/2.4.41 (Ubuntu) Server at colibri-access Port 80</address>
</body></html>"
