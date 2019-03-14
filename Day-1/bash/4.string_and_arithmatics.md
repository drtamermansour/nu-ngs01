<h2 align="center">Loops</h2>

```bash
cd creatures
for filename in basilisk.dat unicorn.dat
do
    head -n 3 $filename
done
```

<h3><b>Can we save the output to file?</b></h3>

```bash
for species in cubane ethane methane
do
    for temperature in 25 30 37 40
    do
       mkdir $species-$temperature
    done
done
```

```bash
while read line;
do
    echo $line | cut -c1-4 ; 
done < basilisk.dat > outputfile.txt
```

```bash
cat basilisk.dat  | while read line;
do
    echo $line | cut -c1-4 ;  
done > outputfile.txt
```

<hr>

<h2 align="center">integer arithmetic expressions</h2>


<b>Bash is restricted to integer arithmetic (ok, there is a way to go around that but wait :)</b>


```bash
k=2; j=5; x=$((k+j)); echo $x
k=2; j=5; let x=k+j; echo $x
declare -i k=2; declare -i j=7; declare -i x; x=k+j; echo $x
```

<h3><b>Useful in conditions for arithmetic conditions</b></h3>

```bash
if (( x == 0 )); …
while (( z > 1000 )); …
if (( p >= 0 && p <= 100 )); …
```

<h3><b>Useful for loop increment/decrement</b></h3>

``bash
let i++
let i--
```

<h3><b>bc for arithmetic expressions</b></h3>

```bash
echo '6.5 / 2.7' | bc
echo 'scale=3; 6.5/2.7' | b
```

```bash
#!/bin/bash
echo "scale=3;
var1 = 6.5 / 2.7;
var2 = 14 * var1;
var2 " \
| bc
```

<hr>

<h2 align="center">String operations</h2>

```bash
var="Welcome to the Bash Scripting"

#Identify String Length inside Bash Shell Script
echo ${#var}
#Extract substring from $string at $position
echo ${var:15}
echo ${var:15:4}

#Shortest Substring Match - chop off prefixes
echo ${var#*t}
echo ${var#Welcome to the }
#Shortest Substring Match - chop off suffixes.
echo ${var%t*}
echo ${var% Scripting}

var="Welcome to the Bash Scripting"

#Longest Substring Match - chop off prefixes
echo ${var##*t}
#Longest Substring Match - chop off suffixes.
echo ${var%%t*}

#Find and Replace String Values inside Bash Shell Script
echo ${var/the Bash/shell}
#Replace all the matches
echo ${var//t/T}
```

<h2 align="center">Shell Scripts</h2>

<ul>
<li><b>You can write the commands directly or in script</b></li>
<li><b>Comments in Bash scripts start with a #</b></li>
<li><b>Bash scripts start with #!/bin/bash (“shebang”)</b></li>
<li><b>Line breaks are bridged with “\”</b></li>
<li><b>The $ is used to label the variable names</b></li>
<li><b>$1, $2, $3,.. are the arguments passed to a script</b></li>
<li><b>The $$ is the process id (PID)</b></li>
</ul>

<h4><b>Do NOT use Microsoft-like word software as a text editors</b></h4>

<hr><hr><hr>

<h2Make your first script></h2>
<br>
<h2><b>Task:</b>Select lines 11-15 from octane.pdb</h2>

```
cd molecules
nano middle.sh
```

```bash
#!/bin/bash 
head -n 15 octane.pdb | tail -n 5
bash middle.sh
```


<h3>Make your script more generic “more usable”
Task: select lines 11-15 from ANY file</h3>

`nano middle.sh`

```bash
#!/bin/bash 
head -n 15 "$1"| tail -n 5
bash middle.sh octane.pdb
bash middle.sh pentane.pdb
```


<h3>More and more generic
Task: select ANY lines from ANY file </h3>

```bash
#!/bin/bash 
head -n "$2" "$1"| tail -n "$3"
```

<h3><b>Can you run this script to select line 8-10 from pentane.pdb?</b></h3>

<h3><b>Can you put this line in a script?

`wc -l *.pdb | sort -n</b></h3>`

```bash
# Sort filenames by their length.
# Usage: bash sorted.sh one_or_more_filenames
wc -l "$@" | sort -n
```

<p align="center"><a href="#"><img src="./assets/8.png"></a></p>


<b>Task: Check for the no of script input parameters</b>

`nano test_input.sh`

```bash
#!/bin/bash 
If (( $# == 2 )); then echo there are 2 parameters;
else echo The number of parameters is $#; fi
```

`bash test_input.sh octane.pdb`

<hr>






