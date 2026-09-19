# Workshop agenda

Written 2026-09-17, from the venue and timing note sent by the SETAC organisers
on the same day. It settles step 4 of `delivery-format-human.md`, the triage of
the material against the time available, and gives each module a stated number
of minutes.

The participant-facing version, which gives the times and the topic blocks
without the triage, is `workshop-agenda-participants.md`. Every timing decision
is stated in full here and summarised there.

## The time available

The organisers' note fixes four points in the day. Instructors reach Room 303
from 09:00, participants from 09:30, and the workshop begins at 10:00. The
opening ceremony of the conference begins at 16:30, and the workshop is asked to
finish in time for participants to attend it, which puts the end of the day at
16:00.

Lunch is served from 12:00 to 13:30 one level below, in bento boxes that can be
taken back to the room, and the full ninety minutes need not be used. Coffee and
tea are in Room 301, next door, from 09:30 to 11:00 and again from 14:00 to
15:00.

That gives 360 minutes between the start and the finish. Taking 45 minutes for
lunch and two 15-minute breaks leaves **285 minutes of teaching**.

`delivery-format-human.md` budgets 380 minutes from a nine-to-five day. The real
day is 95 minutes shorter, so the triage recorded there as the one decision that
cannot be deferred is now the binding constraint on the agenda rather than a
question to return to.

The position of the breaks follows from the coffee station rather than from the
teaching. The morning station closes at 11:00, so the morning break is placed at
10:45 while it is still open. The afternoon station opens at 14:00, so the
afternoon break falls inside the 14:00 to 15:00 window.

## The agenda

| Time | Minutes | Block | Mode |
|---|---|---|---|
| 09:00 | 30 | Room set-up | presenter only |
| 09:30 | 30 | Arrival, coffee in Room 301, installation check desk | hands-on |
| 10:00 | 15 | Welcome, the shape of the day, what a concentration-response model is | opening deck |
| 10:15 | 30 | Module 1, the software stack | demonstration |
| 10:45 | 15 | Break, coffee in Room 301 | |
| 11:00 | 50 | Module 2, fitting a single model | hands-on |
| 11:50 | 40 | Module 3, toxicity estimation and the model set | hands-on |
| 12:30 | 45 | Lunch, one level below | |
| 13:15 | 45 | Module 4, model averaging and multi-model inference | hands-on |
| 14:00 | 30 | Module 5, response data and statistical distributions | demonstration |
| 14:30 | 15 | Break, coffee and snacks in Room 301 | |
| 14:45 | 30 | Module 6, priors and Bayesian inference | demonstration |
| 15:15 | 35 | Module 7, the worked case study, closing with module 8 shown | walk-through |
| 15:50 | 10 | Close, where the material lives, what to read next | |
| 16:00 | | Finish | |
| 16:30 | | Conference opening ceremony | |

Teaching time totals 285 minutes, which is the whole of the budget above.

Three modes are used in the table and are defined here. **Hands-on** means every
participant runs the code as it is introduced. **Demonstration** means the
presenter runs it and participants watch, with the live script available to run
afterwards. **Walk-through** means the page is read and its output discussed
with no code run in the room.

## Time allocated to each module

Modules 2, 3 and 4 are taught in full, hands-on, for 135 minutes between them.
They deliver four of the six outcomes listed in `index.qmd`: fitting a model and
plotting it, deriving NEC, NSEC and ECx, applying model averaging, and
interpreting the weights. Nothing else on the day substitutes for them, and a
participant who runs nothing else should still run these.

Module 2 gets 50 minutes rather than 40 because it is the first time most
participants will have run `bnec()`, and the first hands-on block in any workshop
absorbs the failures that the software setup did not catch.

Modules 5 and 6 are demonstrated rather than run, for 30 minutes each. Module 5
is the slowest module to render, at 12.9 minutes on a four-core machine
(`delivery-format-human.md`, *The time budget*), and its subject is recognising
which distribution suits an endpoint rather than executing a particular call.
Module 6 is the module most likely to need its time adjusted once the survey
reports how many participants have met Bayesian inference before.

Module 7 is a walk-through for 35 minutes. It is the module participants will
copy when they apply the methods to their own data, which is the sixth outcome,
so it earns its place even with no code run in the room.

Module 8 becomes reading. Ten minutes at the end of the module 7 block show what
a grouped comparison looks like and where the module is. It is the most advanced
material in the set and the least likely to be needed by every participant, and
it is the only module whose removal leaves all six outcomes still covered.

Module 1 is cut to 30 minutes against its own length. The software setup is
pre-work, so module 1 explains what R, Stan, `brms` and `bayesnec` each
contribute rather than installing any of them.

## Decisions taken from the survey

The pre-workshop survey in `participant-survey-questions.md` bore on four of the
decisions above, and all four are settled here. None of them changes the total
of 285 teaching minutes, and none of the block times in the agenda above is
altered. What the responses change is the emphasis within three modules and the
staffing of the arrival window.

Responses closed with 18 of the 50 registered participants answering, a response
rate of 36 per cent. The counts below are out of those 18. They are the only
evidence available before the day, and a block that 18 people report needing is
not evidence about the other 32.

The deployed form differs from the design in `participant-survey-questions.md`.
Questions 14, 15, 17 and 18 (operating system, installation rights, point of
failure, and the pasted `check_setup.R` output) are absent, and a single
question asking whether the toolchain is ready to compile and run `bayesnec`
models stands in their place. There is therefore a count of installation
failures and no diagnosis of any of them.

### The time given to module 6

Question 4, comfort with Bayesian statistics, returned 10 responses of "new to
me", 6 of "know the basic idea" and 2 of "have interpreted posteriors before".
No respondent reported working with Bayesian methods routinely. Question 10,
reading a posterior distribution or a credible interval, returned 11 of "not
really" and 7 of "somewhat", and again no respondent reported being comfortable.
Question 12, informative against weakly informative priors, returned 14 of "new
territory".

Module 6 therefore holds at 30 minutes and does not give 10 minutes to module 7.

### The posterior in module 2

The same three questions decide the second half of question 4's purpose, which
was whether the opening block covers what a posterior is before module 2 reaches
one. It does not, because the opening block is 15 minutes and already covers
the welcome, the agenda and the introduction to concentration-response
modelling.

The explanation is given instead at the point in module 2 where the first fit is
plotted and a credible interval appears on the screen, within module 2's
existing 50 minutes. Module 2 was given 50 rather than 40 to absorb the failures
that the software setup did not catch (above), so the reserve is spent on
whichever of the two arises. Where module 2 overruns, module 4 gives back the
5 minutes, for the reason the next section gives.

### The balance within module 3

Question 7, which toxicity values participants report in their own work,
returned 12 responses of "a mix, depending on the study", 2 of "NOEC only", 2 of
"not sure what these are", 1 of "ECx" and 1 of "NEC or NSEC". No single estimate
has a constituency large enough to favour, so module 3 divides its 40 minutes
evenly between NEC, NSEC and ECx.

Question 9, what distinguishes a NEC model from an ECx model, returned 9
responses of "no", 7 of "roughly" and 2 of "yes". Eight respondents chose
"understanding when NEC and when ECx is appropriate" among their priorities for
the day. Module 3's 40 minutes are therefore the minimum rather than a generous
allocation, and the distinction itself is taught rather than assumed.

Two respondents reported not knowing what NOEC, ECx, NEC and NSEC are, and 9 of
18 cannot distinguish a NEC model from an ECx model. The four terms are used
from module 2 onwards. Whether the opening deck names them before module 2 does
is outstanding (RF, 2026-09-19): the slide defining them was cut from the deck
on the ground that the deck should introduce the presenter rather than pre-empt
module 3, and these responses are the case for one sentence naming the four
going back in.

### Module 8 as reading

Questions 6 and 8, the priorities for the day and the intended application,
decide whether module 8 returns to the taught sequence. Question 8 drew 7 free
text answers, of which one describes a grouped comparison: toxicity thresholds
of biodegradable against conventional polymers, and of different particle
shapes, tested for a difference. The others are transcriptomic endpoints (2),
mixtures of low-concentration chemicals with incomplete curves, a pharmaceutical
dose-response curve, a comparison of Bayesian posterior slopes against
frequentist estimates, and one respondent still designing a study.

One respondent of 18 is not the substantial number that would return module 8 to
the sequence. Module 8 stays as reading, shown for 10 minutes at the end of the
module 7 block.

### Module 4 as the reserve

Question 13, the importance of model averaging against selecting a single best
model, returned 7 responses of "not something I currently need", 5 of
"essential", 3 of "useful" and 3 of "not sure what model averaging is".

Module 4 holds at 45 minutes, because model averaging is two of the six outcomes
in `index.qmd` and 5 respondents report needing it. It is named here as the
block to trim first if the day runs late, because it is the only block whose
subject a majority of respondents report not currently needing.

### Staffing of the arrival window

The toolchain question returned 7 responses of "yes", 10 of "not sure" and 1 of
"no". Eleven of 18 respondents have not confirmed that they can compile and run
a model, and 32 registered participants did not answer at all.

Both helpers therefore work the door from 09:30 rather than one, and the
presenter sets up the projector alone. "Not sure" means `check_setup.R` has not
been run, so the reminder sent before the day asks for it to be run and names
the three failures it reports.

### Points for the pre-workshop message

One respondent asked whether the day works in R or in RStudio. Module 1 uses
Positron, so the message states that any editor is fine and that Positron is
what will be on the projector.

Ten respondents chose "confidence installing and running `bayesnec` myself"
among their priorities, equal first with applying the methods to their own data.
Eight chose comparing Bayesian and frequentist results, which is
`vignettes/drc-reference.qmd`. That page is reference rather than taught, so it
is named in the closing block, and it has to be on the approved site by then
(`CLAUDE.md` section 8).

### Experience with R

Question 1 returned 11 responses of "use it occasionally", 6 of "use it
regularly" and 1 of "never used it". Question 2, prior concentration-response
fitting, returned 8 of "a few times", 7 of "never" and 3 of "routinely".
Question 3 returned 9 responses of "none of these" against the list of `drc`,
`bayesnec`, `brms` and Stan; `drc` was the most used, at 5.

This changes no timing, and it is recorded because the hands-on blocks are
written for competent R users. Half the respondents use R occasionally and have
fitted no concentration-response model in any software, so module 2 is the first
time most of the room will have run `bnec()` and the first time many will have
run a model fit of any kind.

## The room and the set-up

Room 303 is classroom style with a table for each attendee and a head table at
the front, and an HDMI connection is available. The seating answers the question
raised with the organisers in August: hands-on work is feasible for 50 people in
this room, and the agenda above assumes it.

One or two helpers are in the room besides the presenter (RF, 2026-09-17), which
settles the staffing question the agenda was drafted against. Modules 2, 3 and 4
therefore stay hands-on. With 50 participants and two helpers the ratio is one
person to 25, and with one helper it is one to 50; at the second ratio a failure
during module 2 may not be reached before module 3 begins, so a helper should
work the room continuously during the hands-on blocks rather than waiting to be
called.

The helpers are given the fits bundle, the live scripts and the three failures
that `check_setup.R` reports, so that a participant is not sent to the presenter
for something the check already names.

The 09:30 to 10:00 arrival window is the only opportunity to catch a broken
installation, and a helper staffs it. One person runs `check_setup.R` from a USB
stick at the door, with the fits bundle on the same stick, while the presenter
sets up the projector.

The fits bundle is a 42.3 MB release asset that downloads in about half a
minute, so the venue network is unlikely to be the constraint for a participant
downloading alone. Several dozen at once is a different matter, and the USB
sticks cover it.

The remaining points from `delivery-format-human.md`, *The room*, are done
before 10:00: the sidebar collapsed, a code chunk checked for legibility from the
back row, and the participants given the address of the approved build rather
than the development one.
