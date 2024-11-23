# REST API

API = application programming interface 

REST stands for representational state transfer. It is a software architectural style that describes a uniform interface between physically separate components, often across the Internet in a Client-Server architecture. 

REST API = RESTful API

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

