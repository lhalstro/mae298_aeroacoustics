%%
%% This is file `aiaa-tc.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% aiaa.dtx  (with options: `class')
%%
%% Welcome to the AIAA LaTeX package.
%%
%% For further details and support, read the Users Manual, aiaa.pdf.
%%
%% Description: LaTeX package for AIAA technical conference papers.
%% Keywords: LaTeX, class, AIAA, BibTeX, bibliographic-style
%% Author: Bil Kleb, Bill Wood, Erich Knausenberger
%% Maintainer: same
%%
%% Run (PDF)LaTeX on `aiaa.dtx' to obtain a Users Manual and code
%% documentation.  Note: to get detailed code documentation you need to
%% comment out the \OnlyDescription line in `aiaa.dtx'.
%%
%% See `aiaa.ins' and the Users Manual for installation instructions.
%%
%% Developer note: synchronize version and date with class code.
%%
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
%%
%%  `aiaa-tc.cls' - A LaTeX class to create AIAA technical conference papers.
%%
%%  The basic approach is to load the LaTeX article class and various
%%  standard LaTeX packages, and then make modifications and extensions to
%%  conform to AIAA's guidelines.
%%
%%  To use, simply select this class file in your LaTeX document.  For example,
%%
%%    \documentclass{aiaa-tc}
%%
%%  Class options:
%%
%%   handcarry - used for hand carrying a paper to a 'loose papers'
%%               conference.  It adds the conference name to the
%%               top of the first page, paper number to the bottom
%%               of each page, and a copyright notice on the bottom
%%               of the first
%%
%%      submit - produces output typically required for journal submission,
%%               that is, 12pt font and 'double spacing'.
%%
%%     handout - Edward Tufte-style format: footnotes and citations become
%%               sidenotes, column width is reduced to give ~60 characters-
%%               per-line to improve readibility, ragged-right used
%%               for improved readibility, and less-obtrusive sectioning
%%               styles.
%%
%%               Note: you will need to move your bibliography section to
%%               the beginning of your document per instructions for the
%%               bibentry LaTeX package.
%%
%%  plus those accepted by the standard article class such as 'draft'.
%%
%%  Class-specific commands:
%%
%%   \thanksibid{} - Reuses \thanks footnote symbol from previous
%%                   author specified by position in \author command.
%%                   For example, `someone\thanksibid{2}' will reuse the
%%                   footnote symbol from the second author for 'someone'
%%                   without producing a repeated footnote at the bottom.
%%
%%   \AIAApapernumber{} - For use with 'handcarry' option.
%%    \AIAAconference{} - Ditto.
%%     \AIAAcopyright{} - Ditto.  Must be filled with one of the
%%                        various copyright types: \AIAAcopyrightA{year},
%%                        \AIAAcopyrightB{year}{holder}, \AIAAcopyrightC,
%%                        \AIAAcopyrightD{year}, or \AIAAcopyrightDAPL{year}.
%%
%%   \AIAAsubmitinfo{} - Footer information used by 'submit' option.
%%
%%  For support, further details, and a documented version of this class,
%%  please read the Users Manual, aiaa.pdf, and the templates that come
%%  with the AIAA LaTeX distribution.
%%
%%  Wish list:
%%
%%  o fine tune hyperref behavior, especially footnote behavior
%%
%%  o an \appendix command to switch section styles, counters, etc.
%%
%% o leave space for copyright notice so that handcarry option does
%%   not alter pagination.
%%
\NeedsTeXFormat{LaTeX2e}
\def\AIAAversion{3.8}
\ProvidesClass{aiaa-tc}[2004/11/23 v3.8 AIAA technical conference paper]
\newif\if@aiaa@handcarry
\newif\if@aiaa@submit
\newif\if@aiaa@handout
\@aiaa@handcarryfalse
\@aiaa@submitfalse
\@aiaa@handoutfalse
\DeclareOption{handcarry}%
 {\@aiaa@handcarrytrue\@aiaa@submitfalse\@aiaa@handoutfalse}
\DeclareOption{submit}%
 {\@aiaa@submittrue\@aiaa@handcarryfalse\@aiaa@handoutfalse}
\DeclareOption{handout}%
 {\@aiaa@handouttrue\@aiaa@handcarryfalse\@aiaa@submitfalse}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}%
                \typeout{Unknown option: passing ``\CurrentOption"
                         option on to the standard LaTeX article class.}}
\ProcessOptions
\if@aiaa@handout
  \LoadClass[11pt]{article}
  \RequirePackage[margin=1in,includemp, letterpaper,%
                  marginparsep=0.375in,
                  marginparwidth=2in,
                  footskip=23pt]{geometry}
  \RequirePackage{paralist}% for compact lists
  \RequirePackage[strict]{chngpage}% for changing page width
  \RequirePackage{bibentry}%  inline citations
\else
 \if@aiaa@submit
  \LoadClass[letterpaper,11pt]{article}
  \RequirePackage{setspace}% control line spacing
 \else
  % \LoadClass[letterpaper,10pt]{article}
  \LoadClass[letterpaper,12pt,
              margin=1in]{article}
 \fi
 \RequirePackage[margin=1in,top=0.75in,headheight=\baselineskip]%
                {geometry}% page layout
 \RequirePackage{overcite}% superscript, sorted, compressed citation commands
\fi
\RequirePackage{fancyhdr}% simple control of headers and footers
\RequirePackage{lastpage}% to get reference for last page number
\RequirePackage{graphicx}% images
\RequirePackage{array}%    extra row height
\newcommand*{\AIAApapernumber}[1]{\def\aiaa@papernumber{#1}}\AIAApapernumber{}
\newcommand*{\AIAAconference}[1]{\def\aiaa@conference{#1}}\AIAAconference{}
\newcommand*{\AIAAprinting}[1]{\def\aiaa@printing{#1}}\AIAAprinting{}
\newcommand*{\AIAAcopyright}[1]{\def\aiaa@copyright{#1}}\AIAAcopyright{}
\newcommand*{\AIAAsubmitinfo}[1]{\def\aiaa@submitinfo{#1}}\AIAAsubmitinfo{}
\setlength{\intextsep}{0.7\baselineskip \@plus 2\p@ \@minus 1\p@}
\def\fps@figure{htbp}
\def\fps@table{htbp}
\setlength\marginparpush{1.6\parsep}
\if@aiaa@handout
 \newcommand{\aiaa@marginfont}{\normalfont\scriptsize}
\fi
\setlength\abovedisplayskip{3pt plus 1pt minus 1pt}
\setlength\belowdisplayskip{3pt plus 1pt minus 1pt}
\setlength\abovecaptionskip{3pt}
\setlength\belowcaptionskip{3pt}
\newlength{\aiaa@overhang}
\setlength{\aiaa@overhang}{\marginparwidth}
\addtolength{\aiaa@overhang}{\marginparsep}
\def\aiaa@raggedright{% resides in latex.ltx
  \let\\\@centercr\@rightskip\z@ plus 0.08\hsize%
  \rightskip\@rightskip
  \leftskip\z@skip
  \parindent\z@}
\newcommand{\AIAAcopyrightA}[1]{Copyright \copyright\ #1 by the
  American Institute of Aeronautics and Astronautics, Inc. All
  rights reserved.}% rev 01/01/2002
\newcommand{\AIAAcopyrightB}[2]{Copyright \copyright\ #1 by
  #2. Published by the American Institute of Aeronautics
  and Astronautics, Inc.\ with permission.}% rev 01/01/2002
\newcommand{\AIAAcopyrightC}{This material is declared a work of the U.S.
  Government and is not subject to copyright protection in the
  United States.}% rev 01/01/2002
\newcommand{\AIAAcopyrightD}[1]{Copyright \copyright\ #1 by the American
  Institute of Aeronautics and Astronautics, Inc.  The U.S. Government
  has a royalty-free license to exercise all rights under the copyright
  claimed herein for Governmental purposes. All other rights are
  reserved by the copyright owner.}% rev 01/01/2002
\newcommand{\AIAAcopyrightDAPL}[1]{Copyright \copyright\ #1 by the American
  Institute of Aeronautics and Astronautics, Inc.  Under the copyright
  claimed herein, the U.S. Government has a royalty-free license to
  exercise all rights for Governmental purposes. All other rights are
  reserved by the copyright owner.}% rev 01/01/2002
\newcommand\aiaa@makecopyright{%
  \begingroup
   \renewcommand\thefootnote{}
   \if@aiaa@handcarry\footnotetext{\aiaa@copyright}\fi
   \if@aiaa@handout\footnotetext{\scriptsize\aiaa@copyright}\fi
  \endgroup
}
\if@aiaa@handout
 \renewcommand{\thanks}[2][0pt]%
   {\footnotemark%
     \renewcommand\thefootnote{\@fnsymbol\c@footnote}
     \marginpar{\aiaa@marginfont
                \vspace*{#1}
                \textsuperscript{\thefootnote}%
                \hspace{1pt}%
                \aiaa@raggedright #2}%
 }
\fi
\newcounter{aiaa@cthanks}
\newcounter{aiaa@cfootnote}
\newcommand{\thanksibid}[1]{%
  \begingroup
    \setcounter{aiaa@cfootnote}{\value{footnote}}% store original
    \setcounter{aiaa@cthanks}{#1}%
    \renewcommand\thefootnote{\fnsymbol{aiaa@cthanks}}%
    \footnotemark%
    \setcounter{footnote}{\value{aiaa@cfootnote}}% reset to original
  \endgroup
}
\if@aiaa@handout
 \renewcommand\maketitle{\par
  \renewcommand{\and}{\vskip 4pt}
  \begingroup
    \parindent 0pt
    \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
    \def\@makefnmark{\rlap{\@textsuperscript{\normalfont\@thefnmark}}}%
    \long\def\@makefntext##1{\parindent 1em\noindent
            \hb@xt@1.8em{%
                \hss\@textsuperscript{\normalfont\@thefnmark}}##1}%
    \newpage
    \global\@topnum\z@   % Prevents figures from going at top of page.
    \@maketitle
    \thispagestyle{plain}\@thanks
  \endgroup
  \setcounter{footnote}{0}%
  \aiaa@makecopyright
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
 }
\else
 \renewcommand\maketitle{\par
  \begingroup
    \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
    \def\@makefnmark{\rlap{\@textsuperscript{\normalfont\@thefnmark}}}%
    \long\def\@makefntext##1{\parindent 1em\noindent
            \hb@xt@1.8em{%
                \hss\@textsuperscript{\normalfont\@thefnmark}}##1}%
    \newpage
    \global\@topnum\z@   % Prevents figures from going at top of page.
    \@maketitle
    \@thanks
  \endgroup
  \setcounter{footnote}{0}%
  \renewcommand\thefootnote{\@alph\c@footnote}% change to alpha style footnote
  \aiaa@makecopyright
  \global\let\thanks\relax
  \global\let\thanksibid\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
 }
\fi
\if@aiaa@handout
 \def\@maketitle{%
  \newpage
  \null
  {\LARGE\bfseries\slshape\raggedright \@title \par}%
  \vskip 10pt%
  {\large\slshape \@author \par}
  \vskip 6pt%
  {\normalsize\slshape AIAA Paper \aiaa@papernumber\par}%
  \vskip 12pt \@plus 3pt \@minus 1pt
 }
\else
 \def\@maketitle{%
  \newpage
  \null
  \begin{center}%
  \if@aiaa@handcarry
    \vspace{-0.25in}
    {\hfill\aiaa@papernumber\par}
    {\sffamily\scriptsize\aiaa@conference}
    \vskip1em
  \else
    \vskip\baselineskip
    \vskip1em
  \fi
  \let \footnote \thanks
    {\LARGE\bfseries \@title \par}%
    \vskip 1.5em%
    {\large
      \lineskip .5em%
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
  \end{center}%
  \par
  \vskip 1.5ex}
\fi
\setlength\hfuzz{2pt}
\fancyhf{}% clear defaults
\renewcommand{\headrulewidth}{0pt}% turn off default headrule
\if@aiaa@handout
  \fancypagestyle{plain}{%
   \fancyhf{}% clear defaults
   \renewcommand{\headrulewidth}{0pt}% turn off default headrule
   \lhead{\sffamily\scriptsize\aiaa@conference}
   \cfoot{\sffamily\scriptsize\aiaa@printing}
  }
\fi
\cfoot{
 \if@aiaa@handout
   \scriptsize\scshape
 \else
   \footnotesize\normalfont
 \fi
 \thepage\ of \pageref{LastPage}\\
 \rule[.2\baselineskip]{0.5in}{0.2pt}\\
 \if@aiaa@submit
   \aiaa@submitinfo
 \else
   \if@aiaa@handout
     \ifx\aiaa@papernumber\@empty
       \relax
     \else
       AIAA Paper \aiaa@papernumber
     \fi
   \else
     % American Institute of Aeronautics and Astronautics
     University of California, Davis
     \ifx\aiaa@papernumber\@empty
       \relax
     \else
       \if@aiaa@handcarry Paper \aiaa@papernumber\fi
     \fi
   \fi
 \fi
}
\pagestyle{fancy}
\setkeys{Gin}{width=\linewidth,totalheight=\textheight,keepaspectratio}
\setlength{\extrarowheight}{2pt}
\if@aiaa@handout
 \renewenvironment{abstract}%
   {\footnotesize\quotation\sffamily\slshape}%
   {\normalfont\endquotation\normalsize}
\else
 \renewenvironment{abstract}%
   {\small\quotation\bfseries}%
   {\normalfont\endquotation\normalsize}
\fi
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\scriptsize\bfseries #1.~~#2}%
  \ifdim \wd\@tempboxa >\hsize
    \scriptsize\bfseries #1.~~#2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
\setlength\belowcaptionskip{5pt plus 1pt minus 1pt}
\if@aiaa@handout
 \setcounter{secnumdepth}{1}
\else
 \renewcommand\thesection      {\@Roman\c@section}
 \renewcommand\thesubsection   {\thesection.\@Alph\c@subsection}
 \renewcommand\thesubsubsection{\thesubsection.\@arabic\c@subsubsection}
 \renewcommand\theparagraph    {\thesubsubsection.\@alph\c@paragraph}
 \renewcommand\thesubparagraph {\theparagraph.\@roman\c@subparagraph}
 \setcounter{secnumdepth}{3}% no counters \paragraph and below
\fi
\if@aiaa@handout\relax\else
 \def\@seccntformat#1{\csname the#1\endcsname.\quad}
\fi
\if@aiaa@handout
 \renewcommand\section{\@startsection {section}{1}{\z@}%
                                      {-2.2ex \@plus -0.5ex \@minus -0.1ex}%
                                      {1.1ex \@plus.1ex}%
                                      {\normalfont\normalsize%
                                       \Large\bfseries\slshape}}
 \renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                        {-1.625ex\@plus -0.5ex \@minus -0.1ex}%
                                        {0.75ex \@plus 0.1ex}%
                                        {\normalfont\normalsize\large\slshape}}
 \renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                           {-1.5ex\@plus -0.4ex \@minus -0.1ex}%
                                           {0.65ex \@plus 0.1ex}%
                                           {\normalfont\normalsize\slshape}}
 \renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                       {2.25ex \@plus1ex \@minus.2ex}%
                                       {-1em}%
                                       {\normalfont\normalsize\slshape}}
 \renewcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                          {2.25ex \@plus1ex \@minus .2ex}%
                                          {-1em}%
                                          {\normalfont\normalfont}}
\else
 \renewcommand\section{\@startsection {section}{1}{\z@}%
                                      {3.5ex \@plus -1ex \@minus -.2ex}%
                                      {2.3ex \@plus.2ex}%
                                      {\center\normalfont\large\bfseries}}
 \renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                        {-3.25ex\@plus -1ex \@minus -.2ex}%
                                        {1.5ex \@plus .2ex}%
                                        {\normalfont\normalsize\bfseries}}
 \renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                           {-3.25ex\@plus -1ex \@minus -.2ex}%
                                           {1.5ex \@plus .2ex}%
                                           {\normalfont\normalsize\itshape}}
 \renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                       {2.25ex \@plus1ex \@minus.2ex}%
                                       {-1em}%
                                       {\normalfont\normalsize\scshape}}
 \renewcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                          {2.25ex \@plus1ex \@minus .2ex}%
                                          {-1em}%
                                          {\normalfont\normalfont}}
\fi
\if@aiaa@handout\relax\else
 \renewcommand{\@biblabel}[1]{$^{#1}$}
 \renewenvironment{thebibliography}[1]%
     {\section*{\refname\@mkboth{\MakeUppercase\refname}%
                                 {\MakeUppercase\refname}}%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\setlength{\leftmargin}{0pt}%
            \settowidth{\labelwidth}{\@biblabel{#1}}%
            \setlength{\itemindent}{\parindent}%
            \advance\itemindent by \labelwidth%
            \setlength{\labelsep}{0.0em}%
            \setlength{\itemsep}{-\smallskipamount}%
            \@openbib@code%
            \usecounter{enumiv}%
            \let\p@enumiv\@empty%
            \renewcommand\theenumiv{\@arabic\c@enumiv}%
            \footnotesize%
            }
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.\@m}
      {\def\@noitemerr
        {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\fi
\if@aiaa@handout
 \newcommand{\aiaa@captionfont}{\scriptsize\sffamily}
 \newcommand{\aiaa@captionfill}{\aiaa@raggedright}
 \long\def\@makecaption#1#2{%
   \vskip\abovecaptionskip
   \sbox\@tempboxa{\aiaa@captionfont\aiaa@captionfill #1: #2}%
   \ifdim \wd\@tempboxa >\hsize
     {\aiaa@captionfont\aiaa@captionfill #1: #2\par}
   \else
     \global \@minipagefalse
     \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
   \fi
    \vskip\belowcaptionskip}
\fi
\if@aiaa@handout
 \renewenvironment{figure*}[1][htbp]
                {\@float{figure}[#1]
                  \begin{adjustwidth}{}{-\aiaa@overhang}%
                   \begin{minipage}{\linewidth}}
                {  \end{minipage}%
                  \end{adjustwidth}%
                 \end@float}
 \renewenvironment{table*}[1][htbp]
                {\@float{table}[#1]
                  \begin{adjustwidth}{}{-\aiaa@overhang}%
                   \begin{minipage}{\linewidth}}
                {  \end{minipage}%
                  \end{adjustwidth}%
                 \end@float}
\fi
\if@aiaa@handout
 \providecommand\newthought[1]{\vspace{2\parsep}{\noindent\scshape #1}}
\fi
\if@aiaa@handout
 \newenvironment{fullwidth}% uses [chngpage.sty]
  {\begin{adjustwidth}{}{-\aiaa@overhang}}%
  {\end{adjustwidth}}
\fi
\if@aiaa@handout
 \newcommand{\sidenote}[2][0pt]%
   {\footnotemark%
    \marginpar{\aiaa@marginfont
               \vspace*{#1}
               \textsuperscript{\thefootnote}%
               \hspace{1pt}%
               \aiaa@raggedright #2}%
   }
\fi
\if@aiaa@handout
 \let\footnote\sidenote
\fi
\if@aiaa@handout
 \newcommand{\marginfigure}[2][0pt]%
  {\marginpar{\def\@captype{figure}
              \aiaa@marginfont
              \vspace*{#1}
              \aiaa@raggedright #2}%
  }
\fi
\if@aiaa@handout
 \newcommand{\margintable}[2][0pt]%
  {\marginpar{\def\@captype{table}
              \aiaa@marginfont
              \vspace*{#1}
              \aiaa@raggedright #2}%
  }
\fi
\AtBeginDocument{%
 \if@aiaa@handout
  \aiaa@raggedright%
  \setlength\parindent{1em}%
  \@ifundefined{wrapoverhang}{}% else wrapfig package loaded
    {\setlength{\wrapoverhang}{\aiaa@overhang}}
  \@ifundefined{captionnamefont}{}% else ccaption package loaded
    {\captionnamefont{\aiaa@captionfont}
     \captiontitlefont{\aiaa@captionfont}
     \captionstyle{\aiaa@raggedright}
    }
 \fi
 \if@aiaa@submit\onehalfspacing\fi
 \@ifundefined{TPT@hsize}{}{%
  \def\tablenotes{%
    \TPT@hsize
    \list{}%
         {\footnotesize\topsep4pt \partopsep\z@skip
          \itemsep2pt \parsep\z@ \itemindent\z@ \leftmargin1.5em
          \def\makelabel##1{\tnote{##1}\hfil}\labelwidth0.5em\labelsep1pt%
          \relax}}
 }%
 \@ifundefined{nomgroup}{}{%
  \RequirePackage{ifthen}% for typesetting groups within nomenclature
  \setlength{\nomitemsep}{-\parsep}
  \renewcommand{\nomgroup}[1]{%
   \medskip%
   \ifthenelse{\equal{#1}{A}}%
    {\relax}%
    {\ifthenelse{\equal{#1}{G}}%
      {\item []\hspace*{-\leftmargin}%
       \emph{Symbols}\hfill\smallskip}%
      {\ifthenelse{\equal{#1}{B}}%
        {\item []\hspace*{-\leftmargin}%
         \emph{Subscripts}\hfill\smallskip}%
        {\ifthenelse{\equal{#1}{C}}%
          {\item []\hspace*{-\leftmargin}%
           \emph{Conventions}\hfill\smallskip}%
          {\ifthenelse{\equal{#1}{T}}%
            {\item []\hspace*{-\leftmargin}%
             \emph{Superscripts}\hfill\smallskip}{}%
          }%
        }%
      }%
    }%
  }% \nomgroup
 }% @ifundefined
}% \AtBeginDocument
\endinput
%%
%% End of file `aiaa-tc.cls'.
