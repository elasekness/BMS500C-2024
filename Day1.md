## Tutorial 1

**Objective:** Successfully log on to Google Cloud Platform (GCP) Console and become familiar with working in a command-line (Linux) environment.

<br>

## Opening a terminal window

We are working on remote servers or virtual machines (VMs) hosted by Google.
To open a terminal window where we will work, we need to log on to the GCP console with our Wadsworth credentials,
start the VM associated with our Project and `ssh` to those servers.
`ssh` is a secure shell protocol for safely connecting to remote services.
To `ssh` you need an account on the server with a username and password as well as the IP address of the server.

- Navigate to the [GCP console](https://console.cloud.google.com) and log on with your Wadsworth credentials
- Once connected to the console, click the Compute Engine link and start your VM by expanding the three dot icon to the right
- Click the SSH button once the VM has started. This will open a terminal window on your VM.

<br>

To perform any operations in a Linux environment, we need to tell the computer what to do by typing specific commands into our terminal window.
In general, the syntax will be: `Command File`,
where `Command` is the function or operation you want to perform on the` File` or `Directory`
that you specify. Sometimes specifying a command is all we need.

<br>

Examine the contents of your directory with:
 
	ls

 > **`ls`** = list

<br>


## Directory structure and navigation

Directories are folders with specific locations on the server.
To navigate from one directory to another, we must specify the location of the directory or its path.
Directories have a forward slash `/` after their names.  Thus to get to a subdirectory within a directory
you would specify the path by stringing together the directory names, separated by `/`. We'll practice below

<br>

First, determine the name and location of your home directory.

	echo $HOME

> **`echo`** = repeat, **`$HOME`** = variable name for your home directory and its location on the server. <br>
> Try the **`echo`** command with other variable names, such as **`$SHELL`** or **`$PATH`**.

<br>

You can also print your working directory (where you are).
 
	pwd

> **`pwd`** = print working directory

<br>

Make a directory within your home directory called `fastq`, where we will store your fastq files.

	mkdir fastq

> **`mkdir`** = make directory

<br>

Now try making three directories at the same time, named `ssu`, `hsp70`, and `actin`.

<br>

Change (move) to a different directory.

	cd fastq

> **`cd`** = change directory. Use **`../`** or **`..`** to move up one directory (back to your home directory). What does **`cd`** alone do? <br>
> You can think of changing directories as physically moving from one directory to another, which means your point of reference has changed. This will become more evident in the upcoming examples. <br>
> For now, try listing the contents of the directory above yours from your home directory. Then **`cd`** to the directory above yours and list the contents of your home directory.

<br>


## Permissions


Directories and files have specific permissions associated with them or things that you are allowed to do to them.
There are three permissions: 1) the ability to read (r) 2) the ability to write (w) and 3) the ability to execute (a script or program; x).
There are three sets of permissions representing 1) the User (you), 2) the group (multiple users),
and 3) Others (Everyone else in the world).

<br>

Examine the permissions of your 'fastq' directory.

	ls -lh fastq

> **`ls -l`** = long list, providing you information on when the **`fastq`** directory was created and its associated permissions. <br>
> the **`-h`** specifies that the output be in human-readable format.

<br>

Change the permissions associated with your `fastq` directory with the 'chmod' command (change mode).

	chmod 775 genomes

> Permissions are represented by three-digit (for user, group, and other) octal numbers.
> Here we are allowing the user and group universal permissions (7 = read, write, and execute) and all others
> the ability to read and write only (5).
> For more information on permissions, see: [Linux permissions](https://www.guru99.com/file-permissions.html#linux_file_ownership)

<br>


## Manipulating files (making, (re)moving, editing, and viewing)

There are many ways to make and view files in a Linux OS (operating system).
We can redirect output from a command that prints its output to your screen (STDOUT) to a file instead,
we can generate files on the fly by opening them in a text editor, or we can copy an existing file.
Similarly, we can view and edit files in a text editor or we can print their contents to the screen with various command-line options.
As a general rule, it's always good to examine some of the contents of your file to ensure you've generated the results you want in the
format you want. Or that you are using the correct file and format for downstream applications.

<br>

Redirect STDOUT to a file.

	ls /usr/bin > programs.txt

> The path **`/usr/bin`** specifies the location where various Bash commands are found. When you type a command, **`/usr/bin`** is one of the locations your computer searches to find and execute the command. <br>
> Was **`/usr/bin`** part of your **`$PATH`**? <br>
> Here we are redirecting the STDOUT from the **`ls`** command to a file named **`programs.txt`**. The **`>`** sign is responsible for the redirection.

<br>

Scroll through the contents of your file.

	more programs.txt

> Scroll through line-by-line with the enter key.  Scroll through page-by-page with the space key. <br>
> Do you notice that the file contains some of the commands you have just used? <br>
> Exit with **`control-c`**.

<br>

Display the first ten lines of your file.

	head programs.txt

> **`head`** displays the first ten lines by default but you can specify the number of lines with a flag. <br>
> For example, **`head -200 programs.txt`** will display the first two hundred lines of your file. In general, most
> commands have additional arguments (or flags) associated with them. <br>
> You can see the different usage statements by typing the command with a **`-help`** or **`--help`** option.

<br>

Display the last ten lines of your file.

	tail programs.txt

<br>

Print the entire contents of your file to your screen.

	cat programs.txt

> **`cat`** = concatenate.  The **`cat`** command can also join the contents of multiple files together.

<br>

Now try copying files from one directory to another.  Here we will copy files required for our future amplicon analysis.
We will copy a fasta file containing reference ssu sequences to our `ssu` directory. We will map our fastq files to these sequences and they will 
also serve as a database for performing BLAST homology searches.  You will generate your own reference/database files for actin and hsp70.
The file are located in different subdirectories within the directory entitled `BMS500-2024`.  Like all of your home directories, `BMS500-2024` is itself 
located in the `/home` directory. We will have to specify the path to these files.  We can specify an absolute path, which is the location
of these files with respect to the root directory (i.e. going through the entire filesystem to get to your file) or a relative path, which is
where these files are located with respect to your current working directory (i.e. where you are when you enter `pwd`).

<br>

Let's try using absolute and relative paths.

	cp /home/BMS500-2024/ssu/ssu_contextual.fasta ssu

> **`cp`** = copy. <br>
> Here we used the absolute path to copy the fasta file `ssu_contextual.fasta` to our `ssu` directory.

<br>

Now `cd` into your `ssu` directory and use a relative path to copy the reference fasta file to where you are. This will simply overwrite the exisiting file.

	cp ../../BMS500-2024/ssu/ssu_contextual.fasta .

> Notice that we had to move up two directories to get to `BMS500-2024`. <br>
> Also notice that we must always specify an end location for our copied files but in this case, we are copying the file to our current location, which is specified with a dot `.`.

<br>

View and edit the contents of a file with a text editor. Let's open our reference fasta file to change the accession name of the first sequence name to the accession followed by the species name and separated by underscores instead of spaces.
We will use a more efficient method of modifying our accession names in the future.

	nano ssu_contextual.fasta

> Nano, emacs, vim, and vi are all text editors. <br>
> You can make an empty file on the fly by typing **`nano`** or **`nano filename`**.  This will open a blank text editor screen. <br>
> Save your changes with **`control-o`**. <br>
> To exit the text editor, use **`control-x`**. <br>
> More information on Nano commands can be found here: 
> [Nano](https://www.howtogeek.com/howto/42980/the-beginners-guide-to-nano-the-linux-command-line-text-editor/)

<br>

Rename a file. You can rename a file with the `mv` (move) command.

	mv ssu_contextual.fasta reference.fasta
  mv reference.fasta ssu_contextual.fasta

> **`mv`** = move. Renaming files with `mv` will overwrite any existing file.  You can also mv a file to a different directory. <br>
> Try it: **`mv ssu_contextual.fasta ../`** <br>
> Can you move the file back to your ssu directory?

<br>

Remove a file. Let's remove files we don't need.

	rm programs.txt

> **`rm`** = remove.  Remember, a removed file cannot be restored. <br>
> Can you remove a file from a different directory without having to change directories? <br>
> How would you remove a directory?

<br>


## Notes on working in a Linux environment


* Avoid creating file names and directories with spaces. Use underscores, dashes, or periods instead to separate multi-part names.
* Spaces need to be escaped in Linux (more on that later).  For example, if you tried to make a directory called “my directory”,
* mkdir would make two directories, `my` and `directory`.
* Use autocomplete for speedier typing and to avoid typos.  Autocomplete will fill in the unique part of a command or file.
For example, if I had only one file in my directory that began with a “b,” I could type `b` and then press the tab key to autocomplete the name of the directory.
* Everything in Linux is case sensitive
* Can’t find a command?
Try `which` to see if the command is in your path – whether the command is in a location that the computer searches for executing commands.
* Hit the up-arrow key to recall a command you entered previously.

<br>

## Databases and obtaining sequences


There are several sequence databases – NCBI, JGI, EMBL - that you might encounter.
You might want to explore each of these to familiarize yourself with the resources they offer.
We will focus on NCBI.  Our goals are to download reference genes for actin and hsp70 and place them in their respective directories.
We can accomplish this several ways.  Here will will download sequences from the webpage and download sequences directly to our VMs with the E-Utilities functions installed on your VM.


Navigate to NCBI’s homepage: [https://www.ncbi.nlm.nih.gov/](https://www.ncbi.nlm.nih.gov/)


> Notice that there are options to submit sequences, download sequences, and even analyze data. <br>
> There is also an option to submit sequences to the BLAST server.
> BLAST is an alignment tool to look for sequences that are similar (a proxy for homology) to your queries, which we will be using later.


Choose 'Nucleotide' under the top left pull-down menu (set to 'All Databases' by default). We can use BOOLEAN search terms (such as AND, OR, and NOT) to make our search more specific.
Type 'actin AND cryptosporidium.'  How many entries are returned? Are they all genes?  Are they all from Cryptospordium (hint: go to the last page of entries for your answer)?.  
We can use the filtering options to the left of the page to make our results even more specific. Limit the results to 'protists' under 'Species' and 'INSDC (GENBANK)' under 'Source databases.'
To ensure we exclude genome assemblies, we can specify a sequence length range. Use the 'custom range' option under 'Sequence length' to specify a sequence length between 400-2000 nucleotides.
We can download these sequences to our computers using the 'send to' option at the top right of the page.  Choose 'Complete record', 'file' as the destination, and 'fasta' as the format.
This would be fine if we were working on our own computers but to use this file on our VMs we have to extra steps.  Transferring this file to the storage bucket associated with the VM and then copying the file from the bucket to the VM.
An easier method to obtain these sequences would be to use one of NCBI's tools for interacting with their databases.

<br>

Return to your VM terminal and type:

	efetch --help

> This brings up a long menu of options for the efetch tool, which can be used to download a variety of data in different formats from NCBI. <br>
> esearch and efetch are part of Entrez Direct suite of utilities ([https://www.ncbi.nlm.nih.gov/books/NBK179288/](https://www.ncbi.nlm.nih.gov/books/NBK179288/)) that allow you to search NCBI databases from your terminal window.
> Relevant arguments are typically the database `-db` we want to search, the format `-format` of the data, and the `-id` of our query, which is the accession of a sequence of interest.

<br>

If we don't have a list of accessions, we will have to search the nucleotide database first and then pass that information to the efetch command, which will download those sequences.
Let's search for Cryptosporidium actin sequences and redirect STDOUT to a file.  First, perform an `esearch` search without piping to `efetch` to see the format of the results.

	esearch -db nuccore -query "actin AND Cryptosporidium [ORGN]"

 > Can you determine how many sequences the esearch command found in the nucleotide database?

Now pass the results of `esearch` to `efetch` to download the sequences.

	esearch -db nuccore -query "actin AND Cryptosporidium [ORGN]" | efetch -format fasta > actin_contextual.fasta

 > Remember that STDOUT is output from a command that gets printed to your screen while the `>` symbol redirects this output to a file.
 > If you look at the definition lines of the sequences or scroll through the file, you'll notice some pretty long sequences.  These are the contigs from various Cryptosporidium genome assemblies
 > that we had exclued previously.  We'll want to do the same here using our BOOLEAN search terms.


	esearch -db nuccore -query "actin AND Cryptosporidium [ORGN] NOT genome" | efetch -format fasta > actin_contextual.fasta

Try using the same commands to download Cryptosporidium hsp70 sequences.

<br>

Although we aren't working with eukaryotic or prokaryotic genomes, it's worth mentioning that there is a command-line
tool to download these as well: the **`datasets`** command-line tool. A help menu will appear if you type **`datasets`** without any arguments.  Typing **`datasets download`** will give you additional information on how to use this command, which shows the option to download a genome by its accession.  


As you can see, there are usually multiple ways to solve a problem in bioinformatics.

<br>


## Manipulating data: Parsing files, modifying content, and piping


Often, we want to extract information from and/or alter the contents of a file. We might also want to determine some basic features of the file. For example, we might want to know how many sequences are in our fasta file without having to open and scroll through a file. You will probably use the following commands most frequently to parse and manipulate files – grep, sed, cut, paste, sort and uniq. These commands can perform simple routines such as search and replace but combined with regular expressions, these tools are incredibly powerful and efficient.

Piping is specified by **`|`** and simply pipes the STDOUT from one command to another so that you can string multiple operations together on one line. Sometimes the most challenging bioinformatics operations are wrangling your data into the proper format.

<br>


`cd` back to your 'actin' directory and count how many sequences are in the 'actin_contextual.fasta' file.

	grep -c ">" actin_contextual.fasta

> **`grep`** = global regular expression print.  Grep searches a file line-by-line for patterns that you specify and returns every line containing that pattern. <br>
> The **`-c`** option counts the number of lines that contain the search pattern instead of returning the lines. <br>
> Try **`grep`** without the **`-c`** argument to see the difference. <br>
> Combined with metacharacters, **`grep`** is a powerful way to search your document for complicated patterns.

<br>

Grab the first 5 header lines from your fasta file with `grep` and a pipe.

	grep  ">" actin_contextual.fasta | head -5

> Here we are using a pipe, designated by **`|`** to capture the output of grep and pass it to another command (**`head`**).
> Piping is a really useful skill to learn for parsing and processing data more efficiently. <br>
> Note that you can string many pipes together, if necessary. As is the case for most operations conducted in Linux, there are multiple ways to do things. <br>
> Use the manual page for grep to find an alternative way to obtain the first five header lines (**`man grep`**).

<br>

Count the number of sequences in the fasta file using a pipe to `wc` instead.

	grep ">" actin_contextual.fasta | wc -l

 > Here we are passing the output of **`grep`** to the word count command, **`wc`**.  The **`-l`** argument specifies that we want to count lines. <br>
 > Again, there is more than one way to achieve the same outcome in Linux.

<br>

Use **`sed`** to rename the definition lines of your fasta file so that only the accessions remain but first `grep` the definition lines so that you can view your changes more easily.


	grep ">" actin_contextual.fasta | sed "s/\.[0-9]*.*//"

> **`sed`** = stream editor.  **`sed`** is essentially a search and replace function. <br>
> Like **`grep`**, we can search for complicated patterns when we use this command with regular expressions. Unlike `grep`, we can replace these complicated patterns with another. <br>
> The syntax for the search and replace command is **`'s/search pattern/replacement pattern/'`** where the 's' stands for substitute. <br>
> In our example, the fasta file contains definition lines that begin with nucleotide accessions, proceeded by metadata separated from the accessions with a spaces. The format of the metadata varies somewhate among the definition lines but note for future iterations of this exercise that all have an acession immediately proceeded by the species name. <br>
> We can use a regular expression to replace all of the patterns displayed in the sequence names without having to search for each pattern individually. Regular expressions use characters with special meaning (metacharacters). <br>
> Like **`grep`**, **`sed`** will search for your pattern line by line. It will make a replacement once (unless you specify otherwise, see the manual page). <br>
> Here we use the regular expression for a **`\.[0-9]*`** to search for a literal period proceeded by any number. The metacharacters: **`.`** and **`*`**. The `.` represents any character one time and the `*` is a greedy character that represents any character any number of times.
> To specify a literal period and not the special character, we need to exit out of the metacharacter with a backlash, giving us `\.`.  Brackets specify a range of numbers or characters.  Here the regular expression `[0-9]*` indicates
> that we are specifying any number from 0-9 (with the brackets) any number of times (with our asterisk). This will find the last part of our accession in every line.  The next expression `.*` represents any character any number of time, which will encompass the rest of the definition line.
> We replace this entire search term with nothing, leaving only the accessions.<br>
Additional metacharacters and their meanings are listed below.


If we didn't `grep` the definitions lines first,  `sed` would print your entire document to STDOUT with the replacements made, making it hard to find the changes, which are only in the definition lines.

* Do you get the same results if you don't include the `.` in your search pattern?
* How would change all spaces to underscores?


<br>


## More piping


Let's try some more complicated parsing of our data using various Bash commands and pipes. 

## Bash for loops


Bash for loops are basically little shell scripts that can be excecuted from the command line (Bash is the command-line language we are using). Like all loops, they allow you to automate iterative processes. For example, instead of opening 200 hundred fasta files and manually changing the definition lines in each, I can run a for loop that will open each fasta file and make the changes that I specify.

<br>

The basic syntax is:

	for FILE in *common_file_ending; do command $FILE; done

> The interpretation of this code is: <br>
> For every file that ends in some common ending (such as .txt or .gz), perform (do) some command on that file until there are no more files on which to operate, whereby “done” will exit us from the loop. <br>
> The $ in front of FILE indicates that $FILE is a variable, a placeholder which is referring to each file that enters the loop, just as x is a variable that represents 2 in the equation x = 2. <br>
> The `for`, `in`, `do`, and `done` are required parts of the for-loop syntax.

<br>

We will try some examples in class.

<br>


## Regular expressions (regex) and special characters (metacharacters)

Regular expressions are search terms that incorporate special characters to make searches more powerful (both broader and more specific).
Metacharacters have a special meaning and include:

- `*` 	Star is a greedy metacharacter meaning match anything any number of times
- `[]` Brackets are often used to specify a range of numbers or letters to include in a search pattern
- `.` 	A period represents any character once
- `?` 	Match one character
- `$`		End of Line
- `^`		Beginning of line

Usually, we escape special characters with a backslash to interpret them literally.

<br>

