# artificial intelligence - AI / chatGPT / etc

xxx move these notes!  leave a link to the course I will likely do next and to the location of my notes

A [Nature intro article](https://www.nature.com/articles/d41586-023-03023-4) by Jeff Leek and others:
- different LLMs are trained on different datasets, and have different "personalities"
- for different uses, we might choose a different LLM

Acronyms:
- "LLM": large language model
- "GPT": Generative Pre-trained Transformer (a type of LLM)

LLM examples include:
- github copilot
- OpenAI chatGPT 3.5 (doesn't include training post-2021)
- OpenAI chatGPT 4   (requires subscription, Maria says it does well for R) Maria pays $20/month. Does hutch have a license? or should we get a group license?
- Bard (from google)
- Claude
- Phind (formerly known as Hello) - designed for software development
- Llama (from Meta = Facebook) 

Strengths and weaknesses
- Bard and Claude might be best for text-y things, not programming things
- Bard will admit when it doesn't know the answer
- chatGPT sometimes MAKES THINGS UP
- Phind provides links to useful sites, but sometimes directly plagiarises

Tips on using LLMs
- be patient
- queries should include what sort of output you want, provide context
- expect some back and forth. Sometimes it takes longer to ask than to just do it.
- test everything - LLMs are fallible. 

Ways to access the LLMs
- copy-paste into the LLM's web interface 
- VScode has a way
- Rstudio does it via the []`chattr` package](https://blogs.rstudio.com/ai/posts/2024-04-04-chat-with-llms-using-chattr/) 


More resources
- DASL course on [using AI for more efficient programming](https://hutchdatascience.org/AI_for_Efficient_Programming/)
- DASL course on [AI considerations applicable to leaders/decision makers](https://hutchdatascience.org/AI_for_Decision_Makers/introduction.html)
- Hutch [AI guidelines](https://centernet.fredhutch.org/u/data-science-lab/ai.html)
- [5 minute video](https://www.youtube.com/watch?v=t7NrkAeosog) on using chatGPT and github copilot within R 
- [Using AI with R](https://rfortherestofus.com/courses/ai) short course form R for the rest of us


Some opinions on when to use AI, and when not to ["AI, LLMs, and BS"](https://datavizf24.classes.andrewheiss.com/resource/ai-bs.html)


[60 min webinar](https://view6.workcast.net/register?cpak=1730348999509149&referrer=OD&utm_campaign=webinars2024&et_rid=35352606&et_cid=5260285) on large language models being used for biological problems 

Maria says github copilot does well for python (older)


Dec 2024:
- Maria/Pravrutha say PerplexityAI is really good for literature search type questions (it gives citations and doesn't make stuff up)
- Maria says Claude Sonnet 3.5 is good for coding

# [Using AI with R](https://rfortherestofus.com/courses/ai) notes

A short course form R for the rest of us

Try to make prompts Simple, Specific, and Use comments

Positron 
- an alternative to Rstudio, newer (still "early stage", in Dec 2024), 
- created by Posit company, the new name for the Rstudio company
- also supports Python (?)
- it is a "fork" of VScode, which is open source


github copilot (owned by Microsoft) in Rstudio
- works within Rstudio (sort of) - code completion works, but chat doesn't.  
- Free account works for most people (2000 code completions and 50 chats per month)
- access from Hutch server Rstudio has been disabled by the admin
- Rstudio - tools - global options, enable copilot, then log in. Probably DO want to index project files (so it can check other files in same project). Probably want automatic suggestions (manual means you need to do a keyboard shortcut)
- start writing code. put in a comment line for what you want copilot to code for you. Tab-complete accepts it.
- code doesn't always work. Can add more comment lines to refine the prompt, delete the original code, and try again
- start writing functions, use tab-complete, AI guesses the rest

`ellmer` package interfaces with a ton of different LLMs (Claude, ChatGPT, google etc etc)
- create API key for each LLM you want to work with
- for chatGPT the login/API key is a different key than you'd use for interacting with chatGPT in the browser. 
- you put the key in an .renviron file
- you use `ellmer::chat_openai()` to start interacting with the chatbot. If you ask it to write code, you then copy-paste output back into a code block

Positron (similar to Rstudio and VScode)
- github copilot doesn't work (MS doesn't want to license it)
- use codeium extension, need to sign up and log in, probably need to create a login token and pass it to Positron
- interact with it similar to using copilot in Rstudio
- codeium might do a bit better than Rstudio
- to update code, highlight code chunk, hit command-I, and instruct codeium to modify code in a particular way

Codeium - an AI 'code accelleration toolkit' that works on a ton of langauges

command-shift-P, type codium, and you get a sidebar that lets you chat with codeium, and it uses context of your code effectively

# Seminar from T (Tara) Templin (University of North Carolina), early 2025

AI can be useful for refining ideas, and for project management.

copilot - good for coding

chatGPT - coding and everything else. Quite good at explaining existing code, for novie coders.

Perplexity - literature searching. You CAN specify that you only want to consider peer-reviewed articles.

ResearchRabbit

Consensus - literature searching

AI models are only as good as the data they are trained on.

Large Language Models (LLMs)
- given the input (question) and the training data, what's the most likely next token?
- Steps, in 'transformer architecture':
  - interpret the input (the user's question)
  - understand the CONTEXT of the question (build "embeddings") - maps tokens from the question onto the training data, figure out what's nearby
  - search for best response
  - craft answer
  - refine answer
  - deliver answer
- autocomplete is a very primitive LLM, similar to AI but with very limited context (e.g. 10 word window)
- 'domain-specific' AIs are usually (at the moment) better than very large LLMs
- in SOME LLMs you can tune the 'temperature' - i.e. the level of creativity or variance
- transformer architecture:
  - can be trained in parallel
  - context can exist between quite distant tokens
  - customizable

Advice on writing questions/prompts for AI:
- ask for creativity, e.g. "give me 20 ideas for ..."
