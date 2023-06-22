sed '/banana/I,+2 d' file #Delete pattern and 2 lines after

sed  's/gr[ae]y/blue/'




perl -0777pe 's/(\n){3,}/\n\n/g' newlines.txt
#Where you should replace \n with whatever newline sequence is appropriate.

#-0777 tells perl to not break each line into its own record, which allows a regex that works across lines to function.

#If you are satisfied with the result, -i causes perl to replace the file in-place rather than output to stdout:

perl -i -0777pe 's/(\n){3,}/\n\n/g' newlines.txt
#You can also do as so: -i~ to create a backup file with the given suffix (~ in this case).

#If slurping the whole file is not acceptable:

perl -ne 'if (/^$/) {$i++}else{$i=0}print if $i<3' newlines.txt
#This prints any line that is not the third (or higher) consecutive empty line. -i works with this the same.

sed -E 's/bla{3,}/BLA/g' #replace with 3 or more patterns



#####

tac | sed '/[A-Za-z\*]{,100}//' | tac #reverses and checks?
tac file | awk '/Query /{c=2} !(c&&c--)' | tac


#####
sed '$!N;/\n.*Query /D;/Query /!P;D' file

#Append the next line (unless the current line is the last line).
#If the appended line contains Query , delete the first line and go again.
#If the first line of the 2 line window contains Query , don't print it.
#Otherwise print the first of the 2 lines, delete it and go again.
#N.B. The appending of the next line is dependant on it not being the last, as the default behaviour of sed is print the pattern space if the N command is called to read passed the end of the file. This allows the last line to treated properly i.e. if the last line contains Query  it will be deleted.



#Replace with sed on a file

#To remove the line and print the output to standard out:
sed '/pattern to match/d' ./infile

#To directly modify the file – does not work with BSD sed:
sed -i '/pattern to match/d' ./infile

#Same, but for BSD sed (Mac OS X and FreeBSD) – does not work with GNU sed:
sed -i '' '/pattern to match/d' ./infile

#To directly modify the file (and create a backup) – works with BSD and GNU sed:
sed -i.bak '/pattern to match/d' ./infile



