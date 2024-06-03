# unix window managers 

## general notes

Choices: tmux, screen, NoMachine

I prefer `screen` - it is a bit of a pain to use, but seems the best way to work around interruptions to the VPN
 
It is worth reading through [this page](https://sciwiki.fredhutch.org/scicomputing/access_methods/#ssh-to-snailfhcrcorg)
 
`NoMachine` is supposed to be a good user-friendly alternative to screen. I tried it a few years ago and did not have good luck, but it sounds like it's what's recommended for most Hutch users.   It was slow for me, and visual resolution was a problem, but maybe it is better now.
 

## screen 

(can put preferences in ~/.screenrc)

screen does not know the aliases that I set up in my login file


1. start webVPN and a single Mac Terminal window
 
2. log in to a rhino node (e.g. jayoung@rhino02.fhcrc.org) - usually you would not bother specifying a particular rhino, and just get on to rhino01 or rhino02 etc. However, if you are using screen, you'll want to specify a particular rhino so that you will be able to restart the same screen session from the same rhino later. 
- I have aliases set up, e.g. `rhino3`
 
3. run the command `screen` to start a screen session. It will be running on whichever rhino you logged in to.
 
4. do whatever you want to do on rhino, e.g. `grabnode` to get onto a gizmo node, run various commands
 
5. if you want a second (third, fourth) command line window within the same screen session, use keystrokes "ctrl-a then c" to get a new window.  To move between windows in a screen session use "ctrl-a then p" (previous) and "ctrl-a then n" (next)
 
If you KNOW you want to pause a session to come back to it later you can 'detach' the screen: "ctrl-a then d"  (but you don't have to - you can still come back to it later even if you don't do this)
 
When you want to come back to your screen session after a break, log in to the same rhino you were on before and run the command `screen -d -r`  (`-d` means detach from whatever client window it was running on before, `-r` means reattach it in the current window).
 
Use the `exit` command to finally terminate the screen process.  You may need to type "exit" several times, e.g. if you have used ssh or grabnode to get on to a gizmo node.
 
### other screen commands: 

- ctrl-a [ = go into copy mode.  There I can scroll up and down within each window. The default number of lines to scroll through is not many but you can increase it as follows:  if needed, create the file `~/.screenrc`. Edit it to include this line, that should tell `screen` to save the most recent 50000 displayed lines:  `defscrollback 50000`

- ctrl-a A to rename a window
- ctrl-a " to get a menu of running windows

- ctrl-a S = split the window horizontally 
- ctrl-a | = split a screen vertically
- ctrl-a Tab = change between splits on the window
- in an empty window/split: To start a new shell in that window, use Ctrl-a c or use Ctrl-a " to choose an existing session to display in that window.
- ctrl-a X = remove a window

- ctrl-a d = detach the session
- screen -ls  = see what screen sessions are running. need to be on the SAME rhino machine as before
- screen -r = reattach a screen session
- screen -d -r = reattach a screen session that's still attached elsewhere (-d means detach it first)
- screen -r PID = reattach a particular screen session


## tmux 

seems a bit less easy than screen, at first glance

https://danielmiessler.com/study/tmux/
https://www.turnkeylinux.org/blog/tmux-screen-alternative

- start webVPN
- start Mac Terminal window
- rhino01 (= ssh jayoung@rhino01.fhcrc.org)
- tmux (tmux is running on rhino01)
   or tmux new -s janetSession1

- grabquarternode (e.g. I get gizmod16)

- ctrl-b d = detach

- tmux -ls = list running sessions
- "tmux a" or "tmux attach" (to reattach to existing session)

- killall tmux = kill any tmux session

## NoMachine

- start webVPN
- start NoMachine
- open sphinx connection.  If a desktop already exists, click on that.
- to leave a sphinx connection open, go to top-right corner to get NoMachine meu up, click Connection-Disconnect and choose "Disconnect" (leaves the session running).  ("Terminate" actually terminates the session)
- getting disconnected because the computer sleeps is also no big deal.

- keyboard shortcuts within the MATE desktop:   
    alt-tab switches between windows
    alt-F6 toggles between windows of the same application
      (can also edit shortcuts)
    shift-control-V to paste

issues:  copy-paste from mac to nomachine window seems very buggy.

https://teams.fhcrc.org/sites/citwiki/SciComp/Pages/Using%20a%20Macintosh%20to%20Access%20SciComp%20Resources.aspx
https://teams.fhcrc.org/sites/citwiki/SciComp/Pages/Connecting%20to%20Session%20Server%20from%20Mac.aspx
https://teams.fhcrc.org/sites/citwiki/SciComp/Pages/Troubleshooting%20NoMachine.aspx


# specific locations

## Seaford (probably old):

router IP 192.168.1.254
my mac via wifi     192.168.1.85
my mac via ethernet 192.168.1.122

