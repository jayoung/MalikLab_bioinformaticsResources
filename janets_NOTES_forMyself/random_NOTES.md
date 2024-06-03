# programming philosophies

['Shaving the yak'](https://www.pipinghotdata.com/posts/2022-11-02-regular-intentional-and-time-boxed-yak-shaving/)


# REST API

API = application programming interface 

REST stands for representational state transfer. It is a software architectural style that describes a uniform interface between physically separate components, often across the Internet in a Client-Server architecture. 

REST API = RESTful API

# Mac tricks / notes

command-. to show hidden files in the Finder (or perhaps command-shift-.)

function-arrow for page up/down (some keyboards have a dedicated button)

if Finder is slow to connect and list files in directories (especially when working remotely), try [this tip](http://osxdaily.com/2015/04/17/fix-slow-folder-populating-cloudkit-macosx/)

apple-shift-3 prints the window

option-copy can copy without going to end of line

shift-return = line break (without paragraph break)

## Excel

To get a line break within a cell in Excel X for the Mac, use `Apple-Option-Enter`

to switch off text to columns:
 one can clear the DATA|Text_to_column settings -- by typing any text character in a cell then activating DATA|text_to_columns, then say choose the Delimited option and clear all of the Delimiters setting boxes.And EXCEL now has a refreshed setting -- thus not requiring to close it and restarting it.

to make an html link in an Excel cell:

row	link
1	=HYPERLINK("http://www.ncbi.nlm.nih.gov/gene/22947", "DUX4")
2	http://www.ncbi.nlm.nih.gov/gene/22947

## Powerpoint

set default font for a presentation - make a text box, get it in the right font. select it and do control-mouseclick. A menu comes up - choose "set as default text box"

when inserting pdfs into a powerpoint, e.g. a doc created using illustrator, do NOT include transparent boxes.  For some reason that makes powerpoint save the image at low resolution, pixellated.
https://answers.microsoft.com/en-us/msoffice/forum/msoffice_powerpoint-mso_mac/inserted-pdfs-become-pixelated-in-powerpoint/bd5e291d-ab01-4771-ab36-f8bfa95da34a
might also be able to adjust how I save out of illustrator to allow me to use transparency: "My fix was to adjust the transparency flattener resolution when outputting the PDF. Setting it to High Resolution resolved my problem." but it seems simpler not to use it

## Endnote

For odd authors, e.g. International Human Genome Sequencing Consortium, make sure author name ends in ", ", and remove from authors term list, and then Endnote won't try to abbreviate it.

## Adobe Illustrator

fonts: some missing fonts in Illustrator are on the computer, but in the wrong place. AdobePiStd seems needed to read R-generated plots.

cp /Library/Application\ Support/Adobe/PDFL/9.0/Fonts/AdobePiStd.otf /Library/Fonts/


quality compression - from large postscript files to less large jpg. Barb's way:

1. make merged pdf in acrobart professional (if multiple images to do)
2. save as jpg (not jpg2000) - output comes as separate files, one per page
3. open in preview


text backgrounds: 

Hi,

Is there a way to set the background color or a text box, so that the
characters are a different color than the text box?

Thanks,
Dan

Yes. Deselect your text box then use your Direct Select or Group Select
tool to select just the box (no text lines should appear), then apply a fill
color. If it's not working check your Appearance palette, make sure it
says Object instead of Type Object to confirm you're editing only the box
object; see if the Fill is showing your new background color and if not clear
Appearance and try again. You may also apply a stroke. You should inset
you text, of course.

inez


## Hutch JAMF service

summer 2017 - Luna set me up on the [Hutch JAMF service](https://centernet.fredhutch.org/cn/u/center-it/help-desk/jamf-pro.html) 

Enroll your mac
Then start up "Self Service" Application
Then install Apple Enterprise Connect (this syncs mac password and hutch password - it is the little key icon on the toolbar)

Can also install Office and Endnote under the centerwide license (no extra cost to center).  Can install Prism but the center pays a per-computer license fee for this.  Later there will be some sort of doc on the website that helps us figure out what is free or not free to the center.


# Zoom tricks

zoom:  try 'view-annotate' to draw marks on other people's shared screen

# html format

HTML tags are not case sensitive: `<P>` means the same as `<p>` (p is the paragraph tag)
opening tag: `<p>`
closing tag: `</p>`

`<a>` for links.  uses href as an attribute in the tag to specify things that happen being the scenes, e.g. 
`<a href="https://www.w3schools.com">This is a link</a>`

img also uses attributes:
`<img src="img_girl.jpg" width="500" height="600">`

can format paragraphs: 
`<p style="color:red">I am a paragraph</p>`

`<pre>` is pre-formatted text (so a fixed pitch font is used)


# rtf - rich text files

## rtf codes for text and background color

These commands work with the header that boxshade uses:

highlight = background color
cf = text color
```
\highlight1\cf0   black text, white background
\highlight1\cf1   white text, white background
\highlight1\cf2   red text, white background
\highlight1\cf3   green text, white background
\highlight1\cf4   blue text, white background
\highlight1\cf5   cyan text, white background
\highlight1\cf6   magenta text, white background
\highlight1\cf7   yellow text, white background
\highlight1\cf8   dark gray text, white background
\highlight1\cf9   light gray text, white background

\highlight8\cf1   white text, grey background
\highlight0\cf1   white text, black background 
```

## rtf headers (from boxshade)
Here's the header:
```
{\rtf1\ansi\deff0
{\fonttbl{\f0\fmodern Courier New;}}
{\info{\author BOXSHADE}}
{\colortbl
\red0\green0\blue0;\red255\green255\blue255;\red255\green0\blue0;\red0\green255\blue0;\red0\green0\blue255;\red0\green255\blue255;\red255\green0\blue255;\red255\green255\blue0;\red128\green128\blue128;\red192\green192\blue192;}
\paperw11880\paperh16820\margl1000\margr500
\margt910\margb910\sectd\cols1\pard\plain
\fs20
```
## rtf paper sizes:

US letter size:
```
\paperw12240
\paperh15840
```
A4 size:
```
\paperw11880
\paperh16820
```

