%%% jhep.bst %%%
%%% History %%%
%%% 2015feb06 sissa medialab %%%

% UT Physics bibliographic style, ver. 1.2.1 Based on:
%
%X     IEEE Transactions bibliography style (29-Jan-88 version)
%X       numeric labels, order-of-reference, IEEE abbreviations,
%X       quotes around article titles, commas separate all fields
%X       except after book titles and before "notes".  Otherwise,
%X       much like the "plain" family, from which this is adapted.
%X
%X       History
%X       9/30/85	(HWT)	Original version, by Howard Trickey.
%X       1/29/88	(OP&HWT) Updated for BibTeX version 0.99a, Oren Patashnik;
%X       THIS `ieeetr' VERSION DOES NOT WORK WITH BIBTEX 0.98i.
%
%  Modifications: 1) added hypertex support and "archive" and "eprint" fields.
%                 2) parentheses around dates, and no "pp." for article entries
%                 3) "publisher, address" instead of "address: publisher"
%                 4) added "report" field for article entries.
%                 5) particle physics-oriented abbreviations, rather than ieee.
%                 6) added "collaboration" field, as per
%                    Jonathan Flynn' suggestion. SPIRES now supports this field.
%                 7) added support for new arXIv identifier format
%                 8) minor corrections (mg 27/2/2012)
%                 9) arXiv links point to arxiv.org and not to xxx.lanl.gov (08/2014)
%
%    Modified by Jacques Distler, 7/97
%           History: ver 1.0  9/96
%                    ver 1.1 10/96 - added "collaboration" field
%                    ver 1.2  7/97 - added a "\providecommand{\href}[2]{#2}"
%                                    to handle case where \href is not defined
%
% HyperTeX Wizardry:
%
% The following are equivalent:
%   archive  = "hep-th"
%   eprint   = "9605023"
% and
%   eprint   = "hep-th/9605023"
% both produce
%    \href{http://arxiv.org/abs/hep-th/9605023}{{\tt hep-th/9605023}}
% in the bibliographic output at the appropriate point. If you are using a
% hypertex macropackage, like hyperref.sty, this will create a link to Los
% Alamos.
%
% The bibtex output produced by SPIRES, while far from perfect, is pretty
% suitable for use with this style. Indeed, this style was designed with
% SPIRES in mind.

ENTRY
  { address
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
    month
    note
    number
    organization
    pages
    publisher
    school
    series
    title
    type
    volume
    year
    doi
    archive
    eprint
    url
    report
    collaboration
    archivePrefix
    primaryClass
  }
  {}
  { label }

INTEGERS { output.state before.all mid.sentence after.quote after.sentence
		after.quoted.block after.block }

FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.quote :=
  #3 'after.sentence :=
  #4 'after.quoted.block :=
  #5 'after.block :=
}

STRINGS { s t }

FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.quote =
	{ " " * write$ }
	{ output.state after.block =
	    { add.period$ write$
	      newline$
	      "\newblock " write$
	    }
	    { output.state before.all =
		'write$
		{ output.state after.quoted.block =
		    { write$
		      newline$
		      "\newblock " write$
		    }
		    { add.period$ " " * write$ }
		  if$
		}
	      if$
	    }
	  if$
	}
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

FUNCTION {output}
{ %duplicate$ write$ 
    duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}

FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

FUNCTION {output.bibitem}
{ newline$
  "\bibitem{" write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {blank.sep}
{ after.quote 'output.state :=
}

FUNCTION {fin.entry}
{ output.state after.quoted.block =
    'skip$
    'add.period$
  if$
  write$
  newline$
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { output.state after.quote =
	{ after.quoted.block 'output.state := }
	{ after.block 'output.state := }
      if$
    }
  if$
}

FUNCTION {new.sentence}
{ output.state after.block =
    'skip$
    { output.state before.all =
	'skip$
	{ after.sentence 'output.state := }
      if$
    }
  if$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

FUNCTION {new.block.checka}
{ empty$
    'skip$
    'new.block
  if$
}

FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
  if$
}

FUNCTION {new.sentence.checka}
{ empty$
    'skip$
    'new.sentence
  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "{\em " swap$ * "}" * }
  if$
}

INTEGERS { nameptr namesleft numnames }

FUNCTION {format.names}
{
        's :=
        #1 'nameptr :=
        s num.names$ 'numnames :=
        numnames 'namesleft :=
        { namesleft #0 > }
        {
                s nameptr "{f.~}{vv~}{ll}{, jj}" format.name$ 't :=
                %pointing to second name and more
                nameptr #1 >
                {
                        numnames #3 >
                        {
                                " {\it et~al.}" *
                                #0 'namesleft :=
                        }
                        {
                                % if more than two names left
                                namesleft #1 >
                                { ", " * t * }
                                {
                                        numnames #2 >
                                        { "," * }
                                        'skip$
                                        if$

                                        t "others" =
                                        { " {\it et~al.}" * }
                                        { " and " * t * }
                                        if$
                                }
                                if$
                        }
                        if$
                }
                't
                if$
                nameptr #1 + 'nameptr :=
                namesleft #1 - 'namesleft :=
        }
        while$
}

FUNCTION {format.authors}
{
        author empty$
        {
                title empty$
                { "" }
                {
                        "{\it " title * "}" *
                }
                if$
        }
        { author format.names }
        if$
}

FUNCTION {format.eprint.newarXivid}
{ eprint #1 #6 substring$ "arXiv:" =
    {"\href{http://arxiv.org/abs/" eprint #7 #10 substring$ * "}" *
        "{{\tt " * eprint #7 #10 substring$ * "}}" *}
    {"\href{http://arxiv.org/abs/" eprint #1 #10 substring$  * "}" *
        "{{\tt arXiv:" * eprint #1 #10 substring$ * "}}" *}
  if$
}

FUNCTION {format.eprint.oldarXivid}
{ archive empty$
    {"\href{http://arxiv.org/abs/" eprint * "}" *
        "{{\tt " * eprint * "}}" *}
    {"\href{http://arxiv.org/abs/" archive *  "/" * eprint * "}" *
        "{{\tt " * archive * "/" * eprint * "}}" *}
  if$
}


FUNCTION {format.eprint}
{ eprint empty$
    { ""}
        { eprint #1 #6 substring$ "arXiv:" =
          eprint #5 #1 substring$ "." =
          or
              { format.eprint.newarXivid }
              { format.eprint.oldarXivid }
        if$
      }
  if$
}


FUNCTION {format.eprint.paren}
{ eprint missing$ { "" } { eprint empty$ { "" } {"[" format.eprint * "]" *}
    			   if$
			  }
  if$
}



FUNCTION {format.report}
{ report empty$
     { ""}
     { report}
     if$
}



FUNCTION {format.editors}
{ editor empty$
    { "" }
    { editor format.names
      editor num.names$ #1 >
	{ ", eds." * }
	{ ", ed." * }
      if$
    }
  if$
}

FUNCTION {format.title}
{ title empty$
    { "" }
    { "{\it " title "t" change.case$ * "}, " * }
  if$
}

FUNCTION {format.title.p}
{ title empty$
    { "" }
    { "``" title "t" change.case$ * ".''" * }
  if$
}

FUNCTION {n.dashify}
{ 't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
	{ t #1 #2 substring$ "--" = not
	    { "--" *
	      t #2 global.max$ substring$ 't :=
	    }
	    {   { t #1 #1 substring$ "-" = }
		{ "-" *
		  t #2 global.max$ substring$ 't :=
		}
	      while$
	    }
	  if$
	}
	{ t #1 #1 substring$ *
	  t #2 global.max$ substring$ 't :=
	}
      if$
    }
  while$
}

FUNCTION {format.date}
{ year empty$
    { month empty$
	{ "" }
	{ "there's a month but no year in " cite$ * warning$
	  month
	}
      if$
    }
    { month empty$
	'year
	{ month ", " * year * }
      if$
    }
  if$
}

FUNCTION {format.date.paren}
{ year empty$
    { month empty$
	{ "" }
	{ "there's a month but no year in " cite$ * warning$
	  month
	}
      if$
    }
    { "(" year * ")" *
%            month empty$
%	{"(" year * ")" *}
%	{"(" month * ", " * year * ")" *}
%      if$
    }
  if$
}

FUNCTION {format.collaboration}
{
    %" test [" collaboration * " Collaboration] t " * write$
    collaboration empty$
    { "" }
    { " [" collaboration * " Collaboration]" * }
  if$
}


FUNCTION {format.btitle}
{ title emphasize
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}

FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { "vol.~" volume *
      series empty$
	'skip$
	{ " of " * series emphasize * }
      if$
      "volume and number" number either.or.check
    }
  if$
}

FUNCTION {format.number.series}
{ volume empty$
    { number empty$
	{ series field.or.null }
	{ output.state mid.sentence =
	    { "no.~" }
	    { "No.~" }
	  if$
	  number *
	  series empty$
	    { "there's a number but no series in " cite$ * warning$ }
	    { " in " * series * }
	  if$
	}
      if$
    }
    { "" }
  if$
}

FUNCTION {format.edition}
{ edition empty$
    { "" }
    { edition "l" change.case$ "~ed." * }
  if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
	{ #1 'multiresult := }
	{ t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}

FUNCTION {format.pages}
{ pages empty$
    { "" }
    { pages multi.page.check
	{ "pp.~" pages n.dashify * }
	{ "p.~" pages * }
      if$
    }
  if$
}

FUNCTION {format.pages.nopp}
{ pages empty$
    { "" }
    { pages multi.page.check
	{ pages n.dashify  }
	{ pages }
      if$
    }
  if$
}


FUNCTION {format.volume}
{ volume empty$
    { "" }
    { "{\bf " volume * "}" * }
  if$
}

FUNCTION {format.number}
{ number empty$
    { "" }
    { "no.~" number * }
  if$
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages
    { type empty$
	{ "ch.~" chapter * }
	{ type "l" change.case$ chapter tie.or.space.connect }
      if$
      pages empty$
	'skip$
	{ ", " * format.pages * }
      if$
    }
  if$
}

FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { "in " booktitle emphasize *
      editor empty$
	'skip$
	{ " (" * format.editors * ")" * }
      if$
    }
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      output.state after.block =
	{ type "t" change.case$ }
	{ type "l" change.case$ }
      if$
    }
  if$
}

FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
  month empty$ year empty$ note empty$
  and and and and and
    { "all relevant fields are empty in " cite$ * warning$ }
    'skip$
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { "Tech. Rep." }
    'type
  if$
  number empty$
    { "l" change.case$ }
    { number tie.or.space.connect }
  if$
}

FUNCTION {format.paddress}
{ address empty$
    { "" }
    { "(" address * ")" * }
  if$
}

FUNCTION {format.article.crossref}
{ key empty$
    { journal empty$
	{ "need key or journal for " cite$ * " to crossref " * crossref *
	  warning$
	  ""
	}
	{ "in { " journal * "\/}" * }
      if$
    }
    { "in " key * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {format.crossref.editor}
{ editor #1 "{vv~}{ll}" format.name$
  editor num.names$ duplicate$
  #2 >
    { pop$ " et~al." * }
    { #2 <
	'skip$
	{ editor #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
	    { " et~al." * }
	    { " and " * editor #2 "{vv~}{ll}" format.name$ * }
	  if$
	}
      if$
    }
  if$
}

FUNCTION {format.book.crossref}
{ volume empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      "In "
    }
    { "Vol.~" volume *
      " of " *
    }
  if$
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
	{ series empty$
	    { "need editor, key, or series for " cite$ * " to crossref " *
	      crossref * warning$
	      "" *
	    }
	    { "{\em " * series * "\/}" * }
	  if$
	}
	{ key * }
      if$
    }
    { format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {format.incoll.inproc.crossref}
{ editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
	{ booktitle empty$
	    { "need editor, key, or booktitle for " cite$ * " to crossref " *
	      crossref * warning$
	      ""
	    }
	    { "in {\em " booktitle * "\/}" * }
	  if$
	}
	{ "in " key * }
      if$
    }
    { "in " format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {first.page}
{ 't :=
  ""
    {
            t empty$ not
            t #1 #1 substring$ "-" = not
            and
    }
    {
            t #1 #1 substring$ *
            t #2 global.max$ substring$ 't :=
    }
    while$
}

FUNCTION {format.journal.pages}
{
        number empty$
        {
                pages empty$
                { "" }
                { " " pages first.page *}
                if$
        }
        {
                " no." number *
                pages empty$
                { "" * }
                { ", " * pages first.page *}
                if$
        }
        if$
}

FUNCTION {first.doi}
{ 't :=
  ""
    {
            t empty$ not
            t #1 #2 substring$ ", " = not
            and
    }
    {
            t #1 #1 substring$ *
            t #2 global.max$ substring$ 't :=
    }
    while$
}

FUNCTION {second.doi}
{ 't :=
  ""
    { t empty$ not }
    {
            t #1 #2 substring$ ", " =
            {
                    t #3 global.max$ substring$ *
                    "" 't :=
            }
            {
                    t #2 global.max$ substring$ 't :=
            }
            if$
    }
    while$
}

FUNCTION {format.journal.name}
{
        % this gives the volume letter if exists
        volume empty$
        { "" }
        {       %checks if first character is a non digit
                volume #1 #1 substring$ chr.to.int$ #47 >
                volume #1 #1 substring$ chr.to.int$ #58 <
                and
                {
                        journal "JHEP" =
                        {
                                year empty$
                                { " {\bf " }
                                { " {\bf " year #3 #4 substring$ * }
                                if$
                        }
                        { " {\bf " }
                        if$
                        volume * "} " *
                }
                {
                        " " volume #1 #1 substring$ *
                        " {\bf " * volume #2 global.max$ substring$ * "} " *
                }
                if$
        }
        if$
}

FUNCTION {format.journal.year}
{
        year empty$
        { "" }
        { "(" year * ")" * }
        if$
}

FUNCTION {format.journal.all}
{
        journal format.journal.name *
        format.journal.year *
        format.journal.pages *
}

FUNCTION {format.journal}
{
        doi missing$
        {
                format.journal.all
        }
        {
                " \href{http://doi.org/" doi first.doi * "}" *
                "{" * format.journal.all * "}" *
        }
        if$
}

FUNCTION {format.journal.erratum}
{
        note missing$
        { "" }
        {
                %" llll " #21 note text.length$ #5 - + int.to.str$ * write$
                note #2 #7 substring$ "Erratum" =
                {
                        " (Erratum: "
                        doi second.doi empty$
                        {
                                note #11 note text.length$ #11 - substring$ * ")" *
                        }
                        {
                                "\href{http://doi.org/" * doi second.doi * "}{" *
                                note #11 note text.length$ #11 - substring$ * "})" *
                        }
                        if$
                }
                { "" }
                if$
        }
        if$
}

FUNCTION {format.arxiv}
{
        "\href{http://arxiv.org/abs/" eprint * "}" *
        eprint #1 #1 substring$ chr.to.int$ #47 >
        eprint #1 #1 substring$ chr.to.int$ #58 <
        and
        { "{arXiv:" * eprint * " [" * primaryClass * "]}" * }
        { "{[" * eprint * "]}" * }
        if$
}


FUNCTION {article}
{
        output.bibitem
        write$
        format.authors 
        write$
        format.collaboration ", " *
        write$
        eprint empty$
        {
                journal empty$
                { "." }
                {
                        format.journal
                        format.journal.erratum * "." *
                }
                if$
        }
        {
                journal empty$
                { "" }
                {
                        format.journal
                        format.journal.erratum * ", " * 
                }
                if$

                format.arxiv * "." *
        }
        if$
        write$
        %eprint empty$
        %{ "." } 
        %{ format.arxiv * ". " * }
        %if$
        %write$
        newline$
}

FUNCTION {link}
{
        output.bibitem
        write$
        format.authors " (" *
        write$
        format.date ") " *
        write$
        title empty$
        { "." }
        { url empty$
           { ", {\it " title "t" change.case$ * "}." * }
           { ", \href{" url * "}{\it " * title "t" change.case$ * "}." * }
           if$
        }
        if$
        write$
        newline$
}


FUNCTION {phdthesis}
{
        output.bibitem
        write$
        format.authors " (" *
        write$
        format.date "). " *
        write$
        format.title 
        write$
        url empty$
        { "Doctoral dissertation." }
        { "\href{" url * "}{Doctoral dissertation}." * }
        if$
        write$
        newline$
}

FUNCTION {mastersthesis} {phdthesis}

FUNCTION {cpp}
{
        output.bibitem
        write$
        format.authors " (" *
        write$
        format.date "). " *
        write$
        format.title 
        write$
        url empty$
        { "C++ library." }
        { "\href{" url * "}{C++ library}." * }
        if$
        write$
        newline$
}

FUNCTION {lecture}
{
        output.bibitem
        write$
        format.authors " (" *
        write$
        format.date "). " *
        write$
        format.title 
        write$
        url empty$
        { "Lecture series." }
        { "\href{" url * "}{Lecture series}." * }
        if$
        write$
        newline$
}

FUNCTION {private}
{
        output.bibitem
        write$
        format.authors 
        write$
        format.collaboration ", " *
        write$
        "Private communication ("
        write$
        format.date "). " *
        write$
        newline$
}

FUNCTION {default.type} { article }

MACRO {jan} {"Jan."}

MACRO {feb} {"Feb."}

MACRO {mar} {"Mar."}

MACRO {apr} {"Apr."}

MACRO {may} {"May"}

MACRO {jun} {"June"}

MACRO {jul} {"July"}

MACRO {aug} {"Aug."}

MACRO {sep} {"Sept."}

MACRO {oct} {"Oct."}

MACRO {nov} {"Nov."}

MACRO {dec} {"Dec."}

MACRO {nup} {"Nucl. Phys."}

MACRO {cmp} {"Comm. Math. Phys."}

MACRO {prl} {"Phys. Rev. Lett."}

MACRO {pl} {"Phys. Lett."}

MACRO {rmp} {"Rev. Mod. Phys."}

MACRO {ijmp} {"Int. Jour. Mod. Phys."}

MACRO {mpl} {"Mod. Phys. Lett."}

MACRO {pr} {"Phys. Rev."}

MACRO {jhep} {"JHEP"}

READ

STRINGS { longest.label }

INTEGERS { number.label longest.label.width }

FUNCTION {initialize.longest.label}
{ "" 'longest.label :=
  #1 'number.label :=
  #0 'longest.label.width :=
}

FUNCTION {longest.label.pass}
{ number.label int.to.str$ 'label :=
  number.label #1 + 'number.label :=
  label width$ longest.label.width >
    { label 'longest.label :=
      label width$ 'longest.label.width :=
    }
    'skip$
  if$
}

EXECUTE {initialize.longest.label}

ITERATE {longest.label.pass}

FUNCTION {begin.bib}
{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\providecommand{\href}[2]{#2}"
  "\begingroup\raggedright\begin{thebibliography}{" * longest.label  *
  "}" * write$ newline$ }

EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}\endgroup" write$ newline$
}

EXECUTE {end.bib}
