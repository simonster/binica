
MAKE(1) 		FreeBSD General Commands Manual 	       MAKE(1)

NNAAMMEE
     mmaakkee - maintain program dependencies

SSYYNNOOPPSSIISS
     mmaakkee [--BBPPSSeeiikknnqqrrssttvv] [--DD _v_a_r_i_a_b_l_e] [--dd _f_l_a_g_s] [--EE _v_a_r_i_a_b_l_e] [--ff _m_a_k_e_f_i_l_e]
	  [--II _d_i_r_e_c_t_o_r_y] [--jj _m_a_x___j_o_b_s] [--mm _d_i_r_e_c_t_o_r_y] [--VV _v_a_r_i_a_b_l_e]
	  [_v_a_r_i_a_b_l_e=_v_a_l_u_e] [_t_a_r_g_e_t _._._.]

DDEESSCCRRIIPPTTIIOONN
     MMaakkee is a program designed to simplify the maintenance of other programs.
     Its input is a list of specifications describing dependency relationships
     between the generation of files and programs.  The first of `_m_a_k_e_f_i_l_e'
     and `_M_a_k_e_f_i_l_e' that can be found in either the current directory or a
     special object directory (see `_._O_B_J_D_I_R') will be read for this list of
     specifications.  If the file `_._d_e_p_e_n_d' can be found, it is also read (see
     mkdep(1)).

     This manual page is intended as a reference document only.  For a more
     thorough introduction to mmaakkee and makefiles, please refer to _M_a_k_e _- _A
     _T_u_t_o_r_i_a_l.

     The options are as follows:

     --BB      Try to be backwards compatible by executing a single shell per
	     command and by executing the commands to make the sources of a
	     dependency line in sequence.  This is turned on by default unless
	     --jj is used.

     --DD _v_a_r_i_a_b_l_e
	     Define _v_a_r_i_a_b_l_e to be 1, in the global context.

     --dd _f_l_a_g_s
	     Turn on debugging, and specify which portions of mmaakkee are to
	     print debugging information.  Argument _f_l_a_g_s is one or more of
	     the following:

	     _A	     Print all possible debugging information; equivalent to
		     specifying all of the debugging flags.

	     _a	     Print debugging information about archive searching and
		     caching.

	     _c	     Print debugging information about conditional evaluation.

	     _d	     Print debugging information about directory searching and
		     caching.

	     _f	     Print debugging information about the execution of for
		     loops.  Currently a no-op.

	     _g_1      Print the input graph before making anything.

	     _g_2      Print the input graph after making everything, or before
		     exiting on error.

	     _j	     Print debugging information about running multiple
		     shells.

	     _l	     Print commands in Makefiles regardless of whether or not
		     they are prefixed by @ or other "quiet" flags.  Also

		     known as "loud" behavior.

	     _m	     Print debugging information about making targets, includ-
		     ing modification dates.

	     _s	     Print debugging information about suffix-transformation
		     rules.

	     _t	     Print debugging information about target list mainte-
		     nance.

	     _v	     Print debugging information about variable assignment.

     --EE _v_a_r_i_a_b_l_e
	     Specify a variable whose environment value (if any) will override
	     macro assignments within makefiles.

     --ee      Specify that environment values override macro assignments within
	     makefiles for all variables.

     --ff _m_a_k_e_f_i_l_e
	     Specify a makefile to read instead of the default `_m_a_k_e_f_i_l_e' and
	     `_M_a_k_e_f_i_l_e'. If _m_a_k_e_f_i_l_e is `-', standard input is read.  Multiple
	     makefiles may be specified, and are read in the order specified.

     --II _d_i_r_e_c_t_o_r_y
	     Specify a directory in which to search for makefiles and included
	     makefiles.  The system makefile directory (or directories, see
	     the --mm option) is automatically included as part of this list.

     --ii      Ignore non-zero exit of shell commands in the makefile.  Equiva-
	     lent to specifying `-' before each command line in the makefile.

     --jj _m_a_x___j_o_b_s
	     Specify the maximum number of jobs that mmaakkee may have running at
	     any one time.  Turns compatibility mode off, unless the _B flag is
	     also specified.

     --kk      Continue processing after errors are encountered, but only on
	     those targets that do not depend on the target whose creation
	     caused the error.

     --mm _d_i_r_e_c_t_o_r_y
	     Specify a directory in which to search for sys.mk and makefiles
	     included via the <...> style.  Multiple directories can be added
	     to form a search path.  This path will override the default sys-
	     tem include path: _/_u_s_r_/_s_h_a_r_e_/_m_k. Furthermore, the system include
	     path will be appended to the search path used for "..."-style in-
	     clusions (see the --II option).

     --nn      Display the commands that would have been executed, but do not
	     actually execute them.

     --PP      Collate the output of a given job and display it only when the
	     job finishes, instead of mixing the output of parallel jobs to-
	     gether.  This option has no effect unless --jj is used too.

     --qq      Do not execute any commands, but exit 0 if the specified targets
	     are up-to-date and 1, otherwise.

     --rr      Do not use the built-in rules specified in the system makefile.

     --SS      Stop processing when an error is encountered.  Default behaviour.
	     This is needed to negate the --kk option during recursive builds.

     --ss      Do not echo any commands as they are executed.  Equivalent to

	     specifying `@@' before each command line in the makefile.

     --tt      Rather than re-building a target as specified in the makefile,
	     create it or update its modification time to make it appear up-
	     to-date.

     --VV _v_a_r_i_a_b_l_e
	     Print mmaakkee's idea of the value of _v_a_r_i_a_b_l_e, in the global con-
	     text.  Do not build any targets.  Multiple instances of this op-
	     tion may be specified; the variables will be printed one per
	     line, with a blank line for each null or undefined variable.

     --vv      Be extra verbose.	For multi-job makes, this will cause file ban-
	     ners to be generated.

     _v_a_r_i_a_b_l_e=_v_a_l_u_e
	     Set the value of the variable _v_a_r_i_a_b_l_e to _v_a_l_u_e.

     There are seven different types of lines in a makefile: file dependency
     specifications, shell commands, variable assignments, include statements,
     conditional directives, for loops, and comments.

     In general, lines may be continued from one line to the next by ending
     them with a backslash (`\'). The trailing newline character and initial
     whitespace on the following line are compressed into a single space.

FFIILLEE DDEEPPEENNDDEENNCCYY SSPPEECCIIFFIICCAATTIIOONNSS
     Dependency lines consist of one or more targets, an operator, and zero or
     more sources.  This creates a relationship where the targets ``depend''
     on the sources and are usually created from them.	The exact relationship
     between the target and the source is determined by the operator that sep-
     arates them.  The three operators are as follows:

     ::	   A target is considered out-of-date if its modification time is less
	   than those of any of its sources.  Sources for a target accumulate
	   over dependency lines when this operator is used.  The target is
	   removed if mmaakkee is interrupted.

     !!	   Targets are always re-created, but not until all sources have been
	   examined and re-created as necessary.  Sources for a target accumu-
	   late over dependency lines when this operator is used.  The target
	   is removed if mmaakkee is interrupted.

     ::::    If no sources are specified, the target is always re-created.  Oth-
	   erwise, a target is considered out-of-date if any of its sources
	   has been modified more recently than the target.  Sources for a
	   target do not accumulate over dependency lines when this operator
	   is used.  The target will not be removed if mmaakkee is interrupted.

     Targets and sources may contain the shell wildcard expressions `?', `*',
     `[]' and `{}'. The expressions `?', `*' and `[]' may only be used as part
     of the final component of the target or source, and must be used to de-
     scribe existing files.  The expression `{}' need not necessarily be used
     to describe existing files.  Expansion is in directory order, not alpha-
     betically as done in the shell.

SSHHEELLLL CCOOMMMMAANNDDSS
     Each target may have associated with it a series of shell commands, nor-
     mally used to create the target.  Each of the commands in this script
     _m_u_s_t be preceded by a tab.  While any target may appear on a dependency
     line, only one of these dependencies may be followed by a creation
     script, unless the `::::' operator is used.

     If the first or first two characters of the command line are `@@' and/or
     `--', the command is treated specially.  A `@@' causes the command not to
     be echoed before it is executed.  A `--' causes any non-zero exit status
     of the command line to be ignored.

VVAARRIIAABBLLEE AASSSSIIGGNNMMEENNTTSS
     Variables in mmaakkee are much like variables in the shell, and, by tradi-
     tion, consist of all upper-case letters.  The five operators that can be
     used to assign values to variables are as follows:

     ==	     Assign the value to the variable.	Any previous value is overrid-
	     den.

     ++==      Append the value to the current value of the variable.

     ??==      Assign the value to the variable if it is not already defined.

     ::==      Assign with expansion, i.e. expand the value before assigning it
	     to the variable.  Normally, expansion is not done until the vari-
	     able is referenced.

     !!==      Expand the value and pass it to the shell for execution and as-
	     sign the result to the variable.  Any newlines in the result are
	     replaced with spaces.

     Any whitespace before the assigned _v_a_l_u_e is removed; if the value is be-
     ing appended, a single space is inserted between the previous contents of
     the variable and the appended value.

     Variables are expanded by surrounding the variable name with either curly
     braces (`{}') or parentheses (`()') and preceding it with a dollar sign
     (`$'). If the variable name contains only a single letter, the surround-
     ing braces or parentheses are not required.  This shorter form is not
     recommended.

     Variable substitution occurs at two distinct times, depending on where
     the variable is being used.  Variables in dependency lines are expanded
     as the line is read.  Variables in shell commands are expanded when the
     shell command is executed.

     The four different classes of variables (in order of increasing prece-
     dence) are:

     Environment variables
	     Variables defined as part of mmaakkee's environment.

     Global variables
	     Variables defined in the makefile or in included makefiles.

     Command line variables
	     Variables defined as part of the command line.

     Local variables
	     Variables that are defined specific to a certain target.  The
	     seven local variables are as follows:

	     _._A_L_L_S_R_C   The list of all sources for this target; also known as
		       `_>'.

	     _._A_R_C_H_I_V_E  The name of the archive file; also known as `_!'.

	     _._I_M_P_S_R_C   The name/path of the source from which the target is to
		       be transformed (the ``implied'' source); also known as
		       `_<'.

	     _._M_E_M_B_E_R   The name of the archive member; also known as `_%'.

	     _._O_O_D_A_T_E   The list of sources for this target that were deemed

		       out-of-date; also known as `_?'.

	     _._P_R_E_F_I_X   The file prefix of the file, containing only the file
		       portion, no suffix or preceding directory components;
		       also known as `_*'.

	     _._T_A_R_G_E_T   The name of the target; also known as `_@'.

	     The shorter forms `_@', `_!', `_<', `_%', `_?', `_>', and `_*' are per-
	     mitted for backward compatibility and are not recommended.  The
	     six variables `_@_F', `_@_D', `_<_F', `_<_D', `_*_F', and `_*_D' are permit-
	     ted for compatibility with AT&T System V UNIX makefiles and are
	     not recommended.

	     Four of the local variables may be used in sources on dependency
	     lines because they expand to the proper value for each target on
	     the line.	These variables are `_._T_A_R_G_E_T', `_._P_R_E_F_I_X', `_._A_R_C_H_I_V_E',
	     and `_._M_E_M_B_E_R'.

     In addition, mmaakkee sets or knows about the following internal variables or
     environment variables:

     _$		A single dollar sign `$', i.e.	`$$' expands to a single dol-
		lar sign.

     _M_A_K_E	The name that mmaakkee was executed with (_a_r_g_v[0]).

     _._C_U_R_D_I_R	A path to the directory where mmaakkee was executed.

     _._O_B_J_D_I_R	A path to the directory where the targets are built.  At
		startup, mmaakkee searches for an alternate directory to place
		target files.  It will attempt to change into this special di-
		rectory and will search this directory for makefiles not found
		in the current directory.  The following directories are tried
		in order:

		1.   ${MAKEOBJDIRPREFIX}/`cwd`
		2.   ${MAKEOBJDIR}
		3.   obj.${MACHINE}
		4.   obj
		5.   /usr/obj/`cwd`

		The first directory that mmaakkee successfully changes into is
		used.  If either MAKEOBJDIRPREFIX or MAKEOBJDIR is defined but
		mmaakkee is unable to change into the corresponding directory,
		then the current directory is used without checking the re-
		mainder of the list.  If they are undefined and mmaakkee is unable
		to change into any of the remaining three directories, then
		the current directory is used.

     _._M_A_K_E_F_L_A_G_S
		The environment variable MAKEFLAGS may contain anything that
		may be specified on mmaakkee's command line.  Its contents are
		stored in mmaakkee's _._M_A_K_E_F_L_A_G_S variable.  Anything specified on
		mmaakkee's command line is appended to the _._M_A_K_E_F_L_A_G_S variable
		which is then entered into the environment as MAKEFLAGS for
		all programs which mmaakkee executes.

     _M_F_L_A_G_S	A synonym for _._M_A_K_E_F_L_A_G_S provided for backward compatibility.

     PWD	Alternate path to the current directory.  Supported if built
		with WANT_ENV_PWD defined.  mmaakkee normally sets `_._C_U_R_D_I_R' to
		the canonical path given by getcwd(3).	However, if the envi-
		ronment variable PWD is set and gives a path to the current
		directory, then mmaakkee sets `_._C_U_R_D_I_R' to the value of PWD in-
		stead.	PWD is always set to the value of `_._O_B_J_D_I_R' for all
		programs which mmaakkee executes.

     _._T_A_R_G_E_T_S	List of targets mmaakkee is currently building.

     _._I_N_C_L_U_D_E_S	See ..IINNCCLLUUDDEESS special target.

     _._L_I_B_S	See ..LLIIBBSS special target.

     _M_A_C_H_I_N_E	Name of the machine architecture mmaakkee is running on, obtained
		from the MACHINE environment variable, or through uname(2) if
		not defined.

     _M_A_C_H_I_N_E___A_R_C_H
		Name of the machine architecture mmaakkee was compiled for, de-
		fined at compilation time.

     _V_P_A_T_H	Makefiles may assign a colon-delimited list of directories to
		_V_P_A_T_H. These directories will be searched for source files by
		mmaakkee after mmaakkee has finished parsing all input makefiles.

     Variable expansion may be modified to select or modify each word of the
     variable (where a ``word'' is whitespace-delimited sequence of charac-
     ters).  The general format of a variable expansion is as follows:

	   {variable[:modifier[:...]]}

     Each modifier begins with a colon and one of the following special char-
     acters.  The colon may be escaped with a backslash (`\').

     EE		 Replaces each word in the variable with its suffix.

     HH		 Replaces each word in the variable with everything but the
		 last component.

     MM_p_a_t_t_e_r_n	 Select only those words that match the rest of the modifier.
		 The standard shell wildcard characters (`*', `?', and `[]')
		 may be used.  The wildcard characters may be escaped with a
		 backslash (`\').

     NN_p_a_t_t_e_r_n	 This is identical to MM, but selects all words which do not
		 match the rest of the modifier.

     QQ		 Quotes every shell meta-character in the variable, so that it
		 can be passed safely through recursive invocations of mmaakkee.

     RR		 Replaces each word in the variable with everything but its
		 suffix.

     SS/_o_l_d___s_t_r_i_n_g/_n_e_w___s_t_r_i_n_g/[gg]
		 Modify the first occurrence of _o_l_d___s_t_r_i_n_g in each word of the
		 variable's value, replacing it with _n_e_w___s_t_r_i_n_g. If a `g' is
		 appended to the last slash of the pattern, all occurrences in
		 each word are replaced.  If _o_l_d___s_t_r_i_n_g begins with a caret
		 (`^'), _o_l_d___s_t_r_i_n_g is anchored at the beginning of each word.
		 If _o_l_d___s_t_r_i_n_g ends with a dollar sign (`$'), it is anchored
		 at the end of each word.  Inside _n_e_w___s_t_r_i_n_g, an ampersand
		 (`&') is replaced by _o_l_d___s_t_r_i_n_g. Any character may be used as
		 a delimiter for the parts of the modifier string.  The an-
		 choring, ampersand, and delimiter characters may be escaped
		 with a backslash (`\').

		 Variable expansion occurs in the normal fashion inside both
		 _o_l_d___s_t_r_i_n_g and _n_e_w___s_t_r_i_n_g with the single exception that a
		 backslash is used to prevent the expansion of a dollar sign

		 (`$'), not a preceding dollar sign as is usual.

     TT		 Replaces each word in the variable with its last component.

     _o_l_d___s_t_r_i_n_g_=_n_e_w___s_t_r_i_n_g
		 This is the AT&T System V UNIX style variable substitution.
		 It must be the last modifier specified.  If _o_l_d___s_t_r_i_n_g or
		 _n_e_w___s_t_r_i_n_g do not contain the pattern matching character _%
		 then it is assumed that they are anchored at the end of each
		 word, so only suffixes or entire words may be replaced.  Oth-
		 erwise _% is the substring of _o_l_d___s_t_r_i_n_g to be replaced in
		 _n_e_w___s_t_r_i_n_g

DDIIRREECCTTIIVVEESS,, CCOONNDDIITTIIOONNAALLSS,, AANNDD FFOORR LLOOOOPPSS
     Directives, conditionals, and for loops reminiscent of the C programming
     language are provided in mmaakkee. All such structures are identified by a
     line beginning with a single dot (`.') character.	The following direc-
     tives are supported:

     ..iinncclluuddee _<_f_i_l_e_>

     ..iinncclluuddee _"_f_i_l_e_"
	     Include the specified makefile.  Variables between the angle
	     brackets or double quotes are expanded to form the file name.  If
	     angle brackets are used, the included makefile is expected to be
	     in the system makefile directory.	If double quotes are used, the
	     including makefile's directory and any directories specified us-
	     ing the --II option are searched before the system makefile direc-
	     tory.

     ..uunnddeeff _v_a_r_i_a_b_l_e
	     Un-define the specified global variable.  Only global variables
	     may be un-defined.

     ..eerrrroorr _m_e_s_s_a_g_e
	     Terminate processing of the makefile immediately.	The filename
	     of the makefile, the line on which the error was encountered and
	     the specified message are printed to standard output and mmaakkee
	     terminates with exit code 1.  Variables in the message are ex-
	     panded.

     Conditionals are used to determine which parts of the Makefile to pro-
     cess.  They are used similarly to the conditionals supported by the C
     pre-processor.  The following conditionals are supported:

     ..iiff [!]_e_x_p_r_e_s_s_i_o_n [_o_p_e_r_a_t_o_r _e_x_p_r_e_s_s_i_o_n _._._.]
	     Test the value of an expression.

     ..iiffddeeff [!]_v_a_r_i_a_b_l_e [_o_p_e_r_a_t_o_r _v_a_r_i_a_b_l_e _._._.]
	     Test the value of a variable.

     ..iiffnnddeeff [!]_v_a_r_i_a_b_l_e [_o_p_e_r_a_t_o_r _v_a_r_i_a_b_l_e _._._.]
	     Test the value of a variable.

     ..iiffmmaakkee [!]_t_a_r_g_e_t [_o_p_e_r_a_t_o_r _t_a_r_g_e_t _._._.]
	     Test the target being built.

     ..iiffnnmmaakkee [!]_t_a_r_g_e_t [_o_p_e_r_a_t_o_r _t_a_r_g_e_t _._._.]
	     Test the target being built.

     ..eellssee   Reverse the sense of the last conditional.

     ..eelliiff [!]_e_x_p_r_e_s_s_i_o_n [_o_p_e_r_a_t_o_r _e_x_p_r_e_s_s_i_o_n _._._.]
	     A combination of `..eellssee' followed by `..iiff'.

     ..eelliiffddeeff [!]_v_a_r_i_a_b_l_e [_o_p_e_r_a_t_o_r _v_a_r_i_a_b_l_e _._._.]

	     A combination of `..eellssee' followed by `..iiffddeeff'.

     ..eelliiffnnddeeff [!]_v_a_r_i_a_b_l_e [_o_p_e_r_a_t_o_r _v_a_r_i_a_b_l_e _._._.]
	     A combination of `..eellssee' followed by `..iiffnnddeeff'.

     ..eelliiffmmaakkee [!]_t_a_r_g_e_t [_o_p_e_r_a_t_o_r _t_a_r_g_e_t _._._.]
	     A combination of `..eellssee' followed by `..iiffmmaakkee'.

     ..eelliiffnnmmaakkee [!]_t_a_r_g_e_t [_o_p_e_r_a_t_o_r _t_a_r_g_e_t _._._.]
	     A combination of `..eellssee' followed by `..iiffnnmmaakkee'.

     ..eennddiiff  End the body of the conditional.

     The _o_p_e_r_a_t_o_r may be any one of the following:

     ||||     logical OR

     &&&&     Logical AND; of higher precedence than `||||'.

     As in C, mmaakkee will only evaluate a conditional as far as is necessary to
     determine its value.  Parentheses may be used to change the order of
     evaluation.  The boolean operator `!!' may be used to logically negate an
     entire conditional.  It is of higher precedence than `&&&&'.

     The value of _e_x_p_r_e_s_s_i_o_n may be any of the following:

     ddeeffiinneedd	 Takes a variable name as an argument and evaluates to true if
		 the variable has been defined.

     mmaakkee	 Takes a target name as an argument and evaluates to true if
		 the target was specified as part of mmaakkee's command line or
		 was declared the default target (either implicitly or explic-
		 itly, see _._M_A_I_N) before the line containing the conditional.

     eemmppttyy	 Takes a variable, with possible modifiers, and evaluates to
		 true if the expansion of the variable would result in an emp-
		 ty string.

     eexxiissttss	 Takes a file name as an argument and evaluates to true if the
		 file exists.  The file is searched for on the system search
		 path (see _._P_A_T_H).

     ttaarrggeett	 Takes a target name as an argument and evaluates to true if
		 the target has been defined.

     An _e_x_p_r_e_s_s_i_o_n may also be an arithmetic or string comparison.  Variable
     expansion is performed on both sides of the comparison, after which the
     integral values are compared.  A value is interpreted as hexadecimal if
     it is preceded by 0x, otherwise it is decimal; octal numbers are not sup-
     ported.  The standard C relational operators are all supported.  If after
     variable expansion, either the left or right hand side of a `====' or `!!=='
     operator is not an integral value, then string comparison is performed
     between the expanded variables.  If no relational operator is given, it
     is assumed that the expanded variable is being compared against 0.

     When mmaakkee is evaluating one of these conditional expressions, and it en-
     counters a word it doesn't recognize, either the ``make'' or ``defined''
     expression is applied to it, depending on the form of the conditional.
     If the form is `..iiffddeeff' or `..iiffnnddeeff', the ``defined'' expression is ap-
     plied.  Similarly, if the form is `..iiffmmaakkee' or `..iiffnnmmaakkee', the ``make''
     expression is applied.

     If the conditional evaluates to true the parsing of the makefile contin-
     ues as before.  If it evaluates to false, the following lines are
     skipped.  In both cases this continues until a `..eellssee' or `..eennddiiff' is
     found.

     For loops are typically used to apply a set of rules to a list of files.
     The syntax of a for loop is:

     ..ffoorr _v_a_r_i_a_b_l_e iinn _e_x_p_r_e_s_s_i_o_n
	     <make-rules>

     ..eennddffoorr
     After the for _e_x_p_r_e_s_s_i_o_n is evaluated, it is split into words.  The iter-
     ation _v_a_r_i_a_b_l_e is successively set to each word, and substituted in the
     mmaakkee--rruulleess inside the body of the for loop.

CCOOMMMMEENNTTSS
     Comments begin with a hash (`#') character, anywhere but in a shell com-
     mand line, and continue to the end of the line.

SSPPEECCIIAALL SSOOUURRCCEESS
     ..IIGGNNOORREE	 Ignore any errors from the commands associated with this tar-
		 get, exactly as if they all were preceded by a dash (`-').

     ..MMAAKKEE	 Execute the commands associated with this target even if the
		 --nn or --tt options were specified.  Normally used to mark re-
		 cursive mmaakkee's.

     ..NNOOTTMMAAIINN	 Normally mmaakkee selects the first target it encounters as the
		 default target to be built if no target was specified.  This
		 source prevents this target from being selected.

     ..OOPPTTIIOONNAALL	 If a target is marked with this attribute and mmaakkee can't fig-
		 ure out how to create it, it will ignore this fact and assume
		 the file isn't needed or already exists.

     ..PPRREECCIIOOUUSS	 When mmaakkee is interrupted, it removes any partially made tar-
		 gets.	This source prevents the target from being removed.

     ..SSIILLEENNTT	 Do not echo any of the commands associated with this target,
		 exactly as if they all were preceded by an at sign (`@').

     ..UUSSEE	 Turn the target into mmaakkee's version of a macro.  When the
		 target is used as a source for another target, the other tar-
		 get acquires the commands, sources, and attributes (except
		 for ..UUSSEE) of the source.  If the target already has commands,
		 the ..UUSSEE target's commands are appended to them.

     ..WWAAIITT	 If special ..WWAAIITT source is appears in a dependency line, the
		 sources that precede it are made before the sources that suc-
		 ceed it in the line.  Loops are not being detected and tar-
		 gets that form loops will be silently ignored.

SSPPEECCIIAALL TTAARRGGEETTSS
     Special targets may not be included with other targets, i.e. they must be
     the only target specified.

     ..BBEEGGIINN	 Any command lines attached to this target are executed before
		 anything else is done.

     ..DDEEFFAAUULLTT	 This is sort of a ..UUSSEE rule for any target (that was used on-
		 ly as a source) that mmaakkee can't figure out any other way to
		 create.  Only the shell script is used.  The ..IIMMPPSSRRCC variable
		 of a target that inherits ..DDEEFFAAUULLTT's commands is set to the
		 target's own name.

     ..EENNDD	 Any command lines attached to this target are executed after

		 everything else is done.

     ..IIGGNNOORREE	 Mark each of the sources with the ..IIGGNNOORREE attribute.  If no
		 sources are specified, this is the equivalent of specifying
		 the --ii option.

     ..IINNCCLLUUDDEESS	 A list of suffixes that indicate files that can be included
		 in a source file.  The suffix must have already been declared
		 with ..SSUUFFFFIIXXEESS; any suffix so declared will have the directo-
		 ries on its search path (see ..PPAATTHH) placed in the _._I_N_C_L_U_D_E_S
		 special variable, each preceeded by a --II flag.

     ..IINNTTEERRRRUUPPTT  If mmaakkee is interrupted, the commands for this target will be
		 executed.

     ..LLIIBBSS	 This does for libraries what ..IINNCCLLUUDDEESS does for include
		 files, except that the flag used is --LL.

     ..MMAAIINN	 If no target is specified when mmaakkee is invoked, this target
		 will be built.  This is always set, either explicitly, or im-
		 plicitly when mmaakkee selects the default target, to give the
		 user a way to refer to the default target on the command
		 line.

     ..MMAAKKEEFFLLAAGGSS  This target provides a way to specify flags for mmaakkee when the
		 makefile is used.  The flags are as if typed to the shell,
		 though the --ff option will have no effect.

     ..NNOOTTPPAARRAALLLLEELL
		 Disable parallel mode.

     ..NNOO__PPAARRAALLLLEELL
		 Same as above, for compatibility with other pmake variants.

     ..OORRDDEERR	 The named targets are made in sequence.

     ..PPAATTHH	 The sources are directories which are to be searched for
		 files not found in the current directory.  If no sources are
		 specified, any previously specified directories are deleted.
		 Where possible, use of ..PPAATTHH is preferred over use of the
		 _V_P_A_T_H variable.

     ..PPAATTHH_s_u_f_f_i_x
		 The sources are directories which are to be searched for suf-
		 fixed files not found in the current directory.  MMaakkee first
		 searches the suffixed search path, before reverting to the
		 default path if the file is not found there.  This form is
		 required for ..LLIIBBSS and ..IINNCCLLUUDDEESS to work.

     ..PPHHOONNYY	 Apply the ..PPHHOONNYY attribute to any specified sources.  Targets
		 with this attribute are always considered to be out of date.

     ..PPRREECCIIOOUUSS	 Apply the ..PPRREECCIIOOUUSS attribute to any specified sources.  If
		 no sources are specified, the ..PPRREECCIIOOUUSS attribute is applied
		 to every target in the file.

     ..SSIILLEENNTT	 Apply the ..SSIILLEENNTT attribute to any specified sources.	If no
		 sources are specified, the ..SSIILLEENNTT attribute is applied to
		 every command in the file.

     ..SSUUFFFFIIXXEESS	 Each source specifies a suffix to mmaakkee. If no sources are
		 specified, any previous specified suffices are deleted.

CCOOMMPPAATTIIBBIILLIITTYY
     Older versions of mmaakkee used MAKE instead of MAKEFLAGS. This was removed
     for POSIX compatibility.  The internal variable _M_A_K_E is set to the same
     value as _._M_A_K_E; support for this may be removed in the future.

     Most of the more esoteric features of mmaakkee should probably be avoided for
     greater compatibility.

EENNVVIIRROONNMMEENNTT
     MMaakkee uses the following environment variables, if they exist: MACHINE,
     MAKE, MAKEFLAGS, MAKEOBJDIR, MAKEOBJDIRPREFIX, and PWD.

FFIILLEESS
     .depend			 list of dependencies
     Makefile			 list of dependencies
     makefile			 list of dependencies
     obj			 object directory
     sys.mk			 system makefile
     /usr/share/mk		 system makefile directory
     /usr/share/doc/psd/12.make  PMake tutorial
     /usr/obj			 default MAKEOBJDIRPREFIX directory.

BBUUGGSS
     The determination of _._O_B_J_D_I_R is contorted to the point of absurdity.

     In the presence of several ..MMAAIINN special targets, mmaakkee silently ignores
     all but the first.

     _._T_A_R_G_E_T_S is not set to the default target when mmaakkee is invoked without a
     target name and no ..MMAAIINN special target exists.

     The evaluation of _e_x_p_r_e_s_s_i_o_n in a test is very simple-minded.  Currently,
     the only form that works is `.if ${VAR} op something' For instance, you
     should write tests as `.if ${VAR} = string' not the other way around,
     which doesn't work.

     For loops are expanded before tests, so a fragment such as:

     .for TMACHINE in ${SHARED_ARCHS}
     .if ${TMACHINE} = ${MACHINE}
	  ...
     .endif
     .endfor
     won't work, and should be rewritten the other way around.

SSEEEE AALLSSOO
     mkdep(1)

     _P_M_a_k_e _- _A _T_u_t_o_r_i_a_l.

HHIISSTTOORRYY
     A mmaakkee command appeared in Version 7 AT&T UNIX.

BSD				March 19, 1994				    11
