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

## Decisions left to the survey

The pre-workshop survey in `participant-survey-questions.md` bears on four of
the decisions above. None of them changes the total, so each is an exchange
between blocks.

Question 4, comfort with Bayesian statistics, decides whether module 6 holds at
30 minutes or takes 10 minutes from module 7, and whether the opening block
needs to cover what a posterior is before module 2 rather than leaving it to
module 6.

Question 7, which toxicity values participants report in their own work, decides
the balance within module 3 between NEC, NSEC and ECx.

Questions 6 and 8, the priorities for the day and the intended application,
decide whether module 8 stays as reading. A substantial number of participants
intending to compare curves across groups would return it to the taught
sequence, and the 30 minutes would come from modules 5 and 6.

Questions 16 to 18, the setup check results, decide how many people work the
arrival window and what they are briefed on. Responses are anonymous, so a
participant reporting a failure cannot be contacted beforehand, and the count of
failures is the only signal available before the day. A high count is the case
for putting both helpers on the door at 09:30 rather than one.

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
