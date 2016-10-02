% Class Information
\def\filename{assignment}
\def\fileversion{1.0}
\def\filedate{2013/09/06}
\NeedsTeXFormat{LaTeX2e}[1996/12/01]
\ProvidesClass{assignment}[%
    \filedate\space\fileversion\space LaTeX document class (Finik).]
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
\LoadClass[a4paper]{article}
% Required Packaged
\RequirePackage{ifthen}
\RequirePackage{graphicx}
\RequirePackage{amsfonts}
\RequirePackage{amsmath}
\RequirePackage{etoolbox}
\makeatletter
\patchcmd{\l@section}
  {\hfil}
  {\leaders\hbox{\normalfont$\m@th\mkern \@dotsep mu\hbox{.}\mkern \@dotsep mu$}\hfill}
  {}{}
\makeatother
\RequirePackage{amssymb}
\RequirePackage{dates}
\RequirePackage{hyperref}
\RequirePackage{fancyhdr}
\RequirePackage{etoolbox}
\RequirePackage{hyperref}
\RequirePackage{subfigure}
\RequirePackage{tabularx}
\RequirePackage{multirow}
\RequirePackage{verbatim}
\RequirePackage{latexsym}
\RequirePackage{index}
\RequirePackage{graphicx}
\RequirePackage{listings}
\RequirePackage{xcolor}
\RequirePackage{courier} % for courier font
\RequirePackage[total={6.5in,9.00in},top=1.5in, left=0.9in, includefoot]{geometry}
\lstset{basicstyle=\ttfamily,
  showstringspaces=false,
  commentstyle=\color{red},
  keywordstyle=\color{black}
}
\renewcommand{\maketitle}{
    \begin{titlepage}
		\thispagestyle{empty}    
		\begin{center}
		\vspace*{0.75cm}
		\includegraphics[scale=0.2]{./mit_logo.png}\\
		\vspace*{4.5cm}
		\LARGE\textmd{\textbf{User's Manual}}\\	
		\vspace{1.25in}
		\LARGE\textmd{\textbf{\hmwkTitle}}\\	
		\textbf{\Large \hmwkProblem} \\ 		
%		\vspace{0.5in} \large{by}\\
		\vspace*{0.5cm}
%		\textbf{\large\hmwkAuthor}\\
		\vspace*{4.0cm}
		\includegraphics[scale=1.1]{./kulik_logo.png}\\
		\vspace*{1.25cm}
		\normalsize\hmwkDate\\ % Insert date here if you want it to appear below
		\end{center}
    \end{titlepage}
    %
    \renewcommand{\thefootnote}{\arabic{footnote}}
    \setcounter{footnote}{1}
    %
}
\pagestyle{fancy}
\linespread{1.1} % Line spacing
% Set up the header and footer
\lhead{v1.0} % Top left header
\chead{\texttt{molSimplify}} % Top center head
\rhead{\includegraphics[scale=0.3]{./kulik_logo.png}} % Top right header
\cfoot{} % Bottom center footer
\rfoot{Page\ \thepage} % Bottom right footer
\renewcommand\headrulewidth{0.4pt} % Size of the header rule
\renewcommand\footrulewidth{0.4pt} % Size of the footer rule
\setlength\parindent{0pt} % Removes all indentation from paragraphs



