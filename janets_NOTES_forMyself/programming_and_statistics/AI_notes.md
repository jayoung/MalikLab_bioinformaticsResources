# artificial intelligence - AI / chatGPT / etc

## General notes

This is a rapidly changing field. Older notes might not apply, models and interfaces and pricing improve all the time.

## Acronyms / terminology

- "LLM": large language model
- "GPT": Generative Pre-trained Transformer (a type of LLM)
- "IDE": Integrated Development Environment

AI - artificial intelligence.

Machine learning is a subset of AI. Can train on large or smaller datasets.  Human intervention needed to fix mistakes.

Deep learning is a subset of machine learning. Requires large datasets for training and use neural nets. Deep learning models can learn from their own mistakes.

## List of different LLMs

- Github Copilot - a "pair programmer". Uses OpenAI's Codex (descendent of GPT-3, and itself gets updated)
- [Claude](https://platform.claude.com/docs/en/home) (from Anthropic) (Claude Sonnet)
- OpenAI chatGPT 3.5 (doesn't include training post-2021)
- OpenAI chatGPT 4  
- Perplexity (for summarizing literature)
- Bard (from google)
- Phind (formerly known as Hello) - designed for software development
- Llama (from Meta = Facebook) 


## Ways to access the LLMs

- Use the LLM's own web interface, and copy-paste into it
- VScode has an interface
- Rstudio does it via the [`chattr` package](https://blogs.rstudio.com/ai/posts/2024-04-04-chat-with-llms-using-chattr/) 
- Positron
  - Databot in Positron (from Posit, uses Claude)


## More resources

Hutch resources I might need to come back to:
- Hutch [AI guidelines](https://centernet.fredhutch.org/u/data-science-lab/ai.html)
- DASL course on [AI considerations applicable to leaders/decision makers](https://hutchdatascience.org/AI_for_Decision_Makers/introduction.html)


Yet to read/view (may or may not be useful):
- [Erick Matsen's blog posts](https://matsen.group/agentic.html)
- YouTube video on machine learning by Christian Dallago of NVIDIA
- DASL course on [using AI for more efficient programming](https://hutchdatascience.org/AI_for_Efficient_Programming/)
- some opinions on when to use AI, and when not to ["AI, LLMs, and BS"](https://datavizf24.classes.andrewheiss.com/resource/ai-bs.html)
- [60 min webinar](https://view6.workcast.net/register?cpak=1730348999509149&referrer=OD&utm_campaign=webinars2024&et_rid=35352606&et_cid=5260285) on large language models being used for biological problems 
- [Coursera links](https://www.coursera.org/articles/ai-vs-deep-learning-vs-machine-learning-beginners-guide)
- description of [AI/ML/DL on AWS](https://aws.amazon.com/compare/the-difference-between-machine-learning-and-deep-learning/)
- [Using AI with R video demo](https://posit.co/workflow-demo/ai-powered-data-science-in-positron/)

Already read, might want to come back to:
- [Using AI with R](https://rfortherestofus.com/courses/ai) short course from R for the rest of us. Notes [BELOW](https://github.com/jayoung/MalikLab_bioinformaticsResources/blob/main/janets_NOTES_forMyself/programming_and_statistics/AI_notes.md#using-ai-with-r-notes).

Maybe less useful?
- [LLM Visualization](https://bbycroft.net/llm). Quite technical - perhaps goes beyond what I care about.

# Reading

## A [Nature intro article](https://www.nature.com/articles/d41586-023-03023-4) by Jeff Leek and others:
- different LLMs are trained on different datasets, and have different "personalities"
- for different uses, we might choose a different LLM

LLM examples include:
- github copilot
- OpenAI chatGPT 3.5 (doesn't include training post-2021)
- OpenAI chatGPT 4  
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

# Talks

## Seminar from T (Tara) Templin (University of North Carolina), early 2025

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

## Erick Matsen talk/demo, Hutch, Dec 3 2025

Erick programs in python using the VSCode IDE, doing version control with git.  He does pair programming using `Claude`, through [Claude's command-line interface](https://code.claude.com/docs/en/quickstart) in VScode.

This [quickstart guide to Claude](https://code.claude.com/docs/en/quickstart) seems very helpful.

(there is a thing called Claude Code interface, but Erick does not like it as much)

Each time you 'onboard' Claude to a project it will read all the code in that git repo.  It's good practise to stay in a virtuous cycle of staging CLEAN code (partly because we don't want Claude to learn from bad code)

Good programming practise is also helpful if we want to use AI:
- functions should be SMALL, each performing one task
- code should be split into multiple files, not one giant one. Erick likes about one viewable page per file.
- ideally we have good tests for our functions and/or for the code as a whole

Using Claude:
- talk to it like you would a junior programmer. Use plain English. 
- ask Claude to make planning documents for each task. E.g. "make a .md document with a plan to refactor code from v1 to v2".  Then you can edit that planning document, or give Claude prompts that result in edits. This is an important step.
- the human stays "in the loop" - we should be involved in the coding process at every step.
- Claude can interact with github (meaning github.com not git). `gh` command-line tools allows you to create issues, for example, or to browse the repo.
- Erick thinks Claude is good at coding in Python. Might be less good in R. Definitely works to rewrite R code into python.  python has a package called pytorch that lets you write your own deep learning models 
- tell Claude whether you want to simply accept its changes, or ask for approocal each time

Using any LLM, even in a more conversational way 
- Think of it as a brainstorming partner
- specifically ask it to be critical and not unduly praiseful (that's their tendency)

Erick's coding cycle:
- 1. make a plan
- 2. edit the plan
- 3. file the plan as a github issue
- 4. reset Claude's context (I think meaning ask it to look at the repo again from a fresh start)
- 5. tell Claude to implement the github issue. Be very specific, saying things like "don't get creative" "stop and ask when there's an error", "create a new git branch"
- 6. use git diffs to look at the proposed changes. Maybe fix stuff.
- 7. ask Claude to make a git pull request
- 8. manually review the pull request 

We can use a SANDBOX environment on our computer if we want to isolate Claude's environment. Like a dev-container.  That would prevent it accessing all your ssh keys, for example. I wrote `devcontainer.json` but I don't know why.

"Sub-agents" - these allow us to start a new process with separate context

We might cycle through `make`, running tests, and fixing errors.

Tips for interacting with Claude on the command-line:
- `!` means run an actual command
- `/` I think is how commands to Claude start. 
- without any leading character, I think stuff you type on the command-line is plain text that's part of your conversation with it. E.g. "what does this project do?" or "explain the folder structure" or "here's a bug where users can submit empty forms - fix it"
- `@` refers to a file
- we typical WATCH as it writes code. Can interrupt with `escape` key if it seems to be going in the wrong direction. Pressing `escape` twice REWINDS to coding.
- can send messages to it while it is working.

Some things I wrote down that I'm not totally sure what they mean
- Planning mode 
- exit / claude --resource (? or something). 

`/context` - I think this is how we see how many tokens are being used in our current tasks and how much context is being stored/ sent up and down to the Claude server.

Custom agents - we can create custom agents, e.g. Erick has one called 'clean code reviewer'. You tell it what coding style you like and it'll make an [appropriate .md file](https://github.com/matsen/fake-comp-bio-project/blob/main/.claude/agents/clean-code-reviewer.md) to instruct the sub-agent how to scan your code. You might invoke the subagent using `@clean-code-reviewer` (it is stored in a .md file)

`CLAUDE.md` file - this is where you can set/store your expectations for how Claude interacts with you. Claude may update it as you tell it more expectations. You create that file using `/init` command. You might tell it to "fail fast with no fallbacks" and you might tell it about your usual coding style/rules.


`/config` allows you to adjust Claude setup:
- there's an auto-compact setting. Erick turns it off because he wants to keep the entire conversation, not a compacted version.
- how much context is best to keep?  More is not always better, as performance can degrade.  Claude's context of ~200,000 is roughly the size of a pdf of a longish paper. Google AI context might be 1,000,000.

Can paste in a screenshot (e.g. a plot) and ask Claude what it thinks of it.

Erick says Claude doesn't do as well with python NOTEBOOKS. GGood for doing analysis tasks, but cannot recode them very well.

Use Claude to help you follow coding best practices;
- test-driven development
- type annotations (i.e. object types)
- perhaps include interactive visualizations to help see that data is clean and processed appropriately. I wrote down "Alterra plot" here but don't know what I mean. Perhaps he was referring to the [Vega-Altair python package](https://altair-viz.github.io/index.html) (a.k.a `altair`) for plotting graphs - seems very similar to ggplot



Erick's follow-up email:
- Thanks for coming today!
- https://matsen.group/agentic.html is the blog series, BlueSky and LinkedIn posts if you want to share socially
- https://github.com/matsen/fake-comp-bio-project is the template we played with today. It has the [devcontainer](https://github.com/matsen/fake-comp-bio-project/tree/main/.devcontainer), [subagents](https://github.com/matsen/fake-comp-bio-project/tree/main/.devcontainer), and [slash command](https://github.com/matsen/fake-comp-bio-project/tree/main/.claude/commands) I described. It also has a nice [statusline](https://github.com/matsen/fake-comp-bio-project/tree/main/.claude/statusline) that shows the context that isn't default and I didn't call out.
- AI policy at the center: https://centernet.fredhutch.org/u/data-science-lab/data-governance/ai/artificial-intelligence-review.html
- Guidelines for using AI agents: https://sciwiki.fredhutch.org/datascience/ai_coding/
- Learn about git: https://hutchdatascience.org/catalog/#reproducible-research



# Online video resources

## Short course from R for the rest of us - [Using AI with R](https://rfortherestofus.com/courses/ai)

Try to make prompts Simple, Specific. Always Use comments

Positron (similar interfaceto Rstudio and VScode)
- an alternative to Rstudio, newer (still "early stage", in Dec 2024), 
- created by Posit company, the new name for the Rstudio company
- also supports Python (?)
- it is a "fork" of VScode, which is open source
- github copilot doesn't work (MS doesn't want to license it)
- use codeium extension, need to sign up and log in, probably need to create a login token and pass it to Positron
- interact with it similar to using copilot in Rstudio
- codeium might do a bit better than Rstudio
- to update code, highlight code chunk, hit command-I, and instruct codeium to modify code in a particular way
- as of Dec 2025, positron offers nice way to use AI for coding assistance. When you get an error message, it'll come up with two buttons - "explain" (explains and offers suggestions) and "fix" (actually fixes code).  [Short video explainer](https://www.youtube.com/watch?v=obivHPe6R2Q). Looks like it can use Anthropic (Claude)

Github copilot (owned by Microsoft) in Rstudio
- works within Rstudio (sort of) - code completion works, but chat doesn't.  
- Free account works for most people (2000 code completions and 50 chats per month)
- access from Hutch server Rstudio has been disabled by the admin
- Rstudio - tools - global options, enable copilot, then log in. Probably DO want to index project files (so it can check other files in same project). Probably want automatic suggestions (manual means you need to do a keyboard shortcut)
- start writing code. put in a comment line for what you want copilot to code for you. Tab-complete accepts it.
- code doesn't always work. Can add more comment lines to refine the prompt, delete the original code, and try again
- start writing functions, use tab-complete, AI guesses the rest

`ellmer` package interfaces with a ton of different LLMs (Claude, ChatGPT, google etc etc). The course video seems to call it `elmer` not `ellmer`, which I don't understand.
- create API key for each LLM you want to work with
- for chatGPT the login/API key is a different key than you'd use for interacting with chatGPT in the browser. 
- you put the key in an `.Renviron` file
- you use `ellmer::chat_openai()` to start interacting with the chatbot. If you ask it to write code, you then copy-paste output back into a code block
- "system prompts" define the way we interact with the AI (e.g. tell it we like tidyverse over base R/lattice), whereas "user prompts" are the actual questions we ask the AI. System prompt could be, e.g. "You are an expert programmer who preferse using the tidyverse". User prompt could be "What is the easiers way to make a histogram?"

Codeium - an AI 'code accelleration toolkit' that works on a ton of langauges

command-shift-P, type codium, and you get a sidebar that lets you chat with codeium, and it uses context of your code effectively

`chores` package (formerly known as `pal`, defined in the video as `pal`). pal developer recommends Claude. pal helps Rstudio interact with the chatbot. Addins menu lets us select the LLM. There is a pal called `ggpal2` that defines a set of useful predefined system prompts `inst/prompts/ggplot2-replace.md`.  
So you install both libraries (`pal/chores` and `ggpal2`), then you use the addins to choose pal, using ggpal2 as the selected pal. Then you start coding - you write words that describes what you want to do, highlight the words, and go to pal-ggplot2, and it'll fill in some code that does what you want to do.

Let's say you have some working base-R code and you higlight it and run pal-ggplot2 - it'll turn your base R code into tidyverse code.

You can make your own pal using `prompt_new`. It makes a new `.md` file and we just write in plain English what we want it to do. (must be a newline at the end of the file)

ALWAYS run the code and make sure it works. It often doesn't. Must always check it's doing what you want it to.

Part 5 - analyze data with AI: that video mostly about analyzing survey responses. Using AI can help get a sense of top themes, especially in the freeform text type responses. Can also e.g. translate from Spanish. It's far from perfect but for something like pulling out themes it seems reasonsable.



### [R with chatGPT and github copilot](https://www.youtube.com/watch?v=t7NrkAeosog) (a 5 minute youtube)

Nov 12 2023

Covers GitHub Copilot and ChatGPT

GitHub Copilot 
- a "pair programmer" works, within IDE
ChatGPT
- a chatbot, back and forth for a converation

Can use them together within Rstudio (need more recent Rstudio, and to enable copilot). 

`chattr` package allows use of chatGPT
- maintained by chattr, sends extra info to the chatbot (e.g. tells it about some R references)
- allows use of chatGPT but many other LLMs

# Databot in Positron, for R-based data analysis

Databot demo videos from Ted Laderas (Hutch DASL): Using Databot on the NHANES dataset: 
- [part 1](https://www.youtube.com/watch?v=qs2GozYUUOk)
- [part 2](https://www.youtube.com/watch?v=lT2J71Jg_ug)
- [part 3](https://www.youtube.com/watch?v=j0KdDMIgcLY)

Databot in Positron:
- command-shift-P gives us the command pane, and there's an option to open databot. 
- Select which LLM to use (choose Claude Sonnet)
- then we type requests in plain english. 
- It generates code but doesn't run it without our permission.  
- You would probably want to copy the useful code chunks into your own qmd/Rmd document. 
- It can make plots and tables and plain-english summaries of what the tables show.
- After it runs it offers options for next steps, like to save its code/findings into markdown files

I have installed the extension but I haven't set it up to connect to a particular LLM yet (because I don't have accounts with the LLMs). 


# Notes on other people's experience

Maria says github copilot does well for python (older)

chatGPT 4   (requires subscription, Maria says it does well for R) Maria pays $20/month. Does hutch have a license? or should we get a group license?

Dec 2024:
- Maria/Pravrutha say PerplexityAI is really good for literature search type questions (it gives citations and doesn't make stuff up)
- Maria says Claude Sonnet 3.5 is good for coding
