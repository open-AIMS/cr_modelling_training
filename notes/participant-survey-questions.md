# Participant survey question set

The pre-workshop survey for the concentration-response modelling course, SETAC
Australasia, September 2026. Written for entry into Microsoft Forms by hand,
because Apps Script is blocked under AIMS organisational policy. The equivalent
Google Apps Script build is kept at `scripts/build_participant_survey.gs` and
remains valid if a Google route becomes available later.

The question set follows the structure proposed by Positron Assistant, with the
software setup section expanded. Four questions were added and two response
types changed; the reasons are given against each.

## Form mechanics

Five points about Microsoft Forms decide how this is entered.

Subtitles are not shown by default. The help text under each question below goes
in the Subtitle field, which is switched on per question through the ellipsis
menu at the lower right of the question, then Subtitle.

A Choice question becomes multi-answer through the Multiple answers toggle, and
gains a free-text option through the Other option item in the ellipsis menu.

Forms has no cap on how many options a respondent selects, so the limit of two on
question 6 is stated in its subtitle and cannot be enforced.

Participants are external to AIMS, so the form is set to accept responses from
anyone, and responses are therefore anonymous. No name or email is collected, so
a participant who reports a failing setup in question 16 cannot be identified,
and cannot be sent instructions for that failure before the day. Read the
consequences under Handling the responses before settling on this.

Anonymous response also rules out a File upload question, since file upload
requires the respondent to sign in with a work or school account.

Long answer fields are capped, and the output of a failing `check_setup.R` run
includes up to forty lines of compiler output, so question 18 asks for the
summary block rather than the whole transcript.

## Form settings

Title: Concentration-response modelling workshop: participant survey

Description: This survey shapes the pacing and emphasis of the workshop. It takes
about four minutes. The software setup section at the end matters most: a
participant whose C++ toolchain is not working cannot fit models on the day, and
there is no troubleshooting slot in the programme.

Under Settings, set who can respond to Anyone can respond, leave the response
receipt off, and set the thank-you message to: Your response has been recorded.
If your setup check did not report READY, we will contact you before the
workshop.

Do not use branching. Question 17 would branch from question 16, but branching in
Forms breaks quietly when questions are later reordered, and the subtitle on
question 17 covers the same ground.

## Section 1. Background and prior experience

### 1. Experience with R

Question: Which best describes your current experience with R?

Type: Choice, single answer, required

Options: Never used it / Use it occasionally / Use it regularly / Expert

### 2. Prior concentration-response fitting

Question: Have you fitted a concentration-response or dose-response model
before, in any software?

Type: Choice, single answer, required

Options: Never / A few times / Routinely

### 3. Software used

Question: Which of these have you used before? Select all that apply.

Type: Choice, multiple answers, Other option on, required

Subtitle: Use Other to name software outside R, for example ToxRat, CETIS or
GraphPad.

Options: drc / bayesnec / brms / Stan directly, through rstan or cmdstanr / None
of these

Changed from the proposed single-answer form. Participants arrive having used
both `drc` and `brms`, and a single answer records only one of them.

### 4. Comfort with Bayesian statistics

Question: How would you describe your comfort with Bayesian statistics?

Type: Choice, single answer, required

Options: New to me / Know the basic idea / Have interpreted posteriors before /
Fit Bayesian models routinely

### 5. Role

Question: What is your day-to-day role?

Type: Choice, single answer, Other option on, required

Options: Regulatory risk assessor / Consultant / Academic researcher / Lab-based
ecotoxicologist

## Section 2. Goals for the day

### 6. Priorities for the day

Question: What are you most hoping to get from the workshop?

Type: Choice, multiple answers, Other option on, required

Subtitle: Select up to two.

Options: Confidence installing and running bayesnec myself / Understanding when
NEC and when ECx is appropriate / Justifying a model choice to a reviewer or
regulator / Comparing Bayesian and frequentist results / Applying the methods to
my own data

Changed from single answer, with a stated limit of two, so that the answers rank
rather than force a single choice.

### 7. Toxicity values reported in their own work

Question: Do you currently need to report NOEC, ECx, or NEC and NSEC values in
your work?

Type: Choice, single answer, required

Options: NOEC only / ECx / NEC or NSEC / A mix, depending on the study / Not sure
what these are

### 8. Intended application

Question: Is there a specific dataset or problem you are hoping to apply this to
afterwards?

Type: Text, long answer, optional

Subtitle: A sentence is enough. Response type, number of concentrations, and what
you need to estimate.

## Section 3. Concept checks

Section description: These set the starting point for the relevant modules.
Answer as you stand now, before the course covers them. All optional.

Drop this whole section if you would rather run these as live checks on the day.

### 9. NEC against ECx

Question: Do you know what distinguishes a NEC model from an ECx model?

Type: Choice, single answer, optional

Options: Yes / Roughly / No

### 10. Reading a posterior

Question: How comfortable are you reading a posterior distribution or a credible
interval?

Type: Choice, single answer, optional

Options: Very / Somewhat / Not really

### 11. Choice of statistical distribution

Question: Have you had to justify a choice of statistical distribution, such as
binomial, beta or Gaussian, for response data?

Type: Choice, single answer, optional

Options: Yes, routinely / Once or twice / Never

### 12. Priors

Question: Do you have a view on informative against weakly informative priors?

Type: Choice, single answer, optional

Options: Have a view / Aware of the question / New territory

### 13. Model averaging

Question: How important is model averaging, against selecting a single best
model, in your own work?

Type: Choice, single answer, optional

Options: Essential / Useful / Not something I currently need / Not sure what
model averaging is

The fourth option was added. A participant who has never met model averaging
cannot place it on a scale from essential to not needed.

## Section 4. Software setup

Section description: The setup instructions and `check_setup.R` were issued with
the joining email. Run the check on the machine you will bring on the day.

### 14. Operating system

Question: Which operating system is on the machine you will bring?

Type: Choice, single answer, Other option on, required

Subtitle: The toolchain differs by platform: Rtools on Windows, the Xcode command
line tools on macOS.

Options: Windows / macOS, Apple silicon / macOS, Intel / Linux

Added. Rtools and the Xcode command line tools fail in different ways, and the
split across the room decides what the reminder email says.

### 15. Installation rights

Question: Can you install software on that machine without an IT request?

Type: Choice, single answer, required

Subtitle: A managed laptop that refuses a compiler install is the most common
reason a participant cannot fit a model on the day. Raise the request now if the
answer is no.

Options: Yes / No, installation needs IT approval / Not sure

Added. This is the failure that cannot be fixed in the room, so it has to surface
weeks in advance.

### 16. Setup check result

Question: Did `check_setup.R` report READY?

Type: Choice, single answer, required

Options: Yes, READY / READY, with warnings / NOT READY / Attempted, but the
script would not run / Have not tried yet

The five options match what `vignettes/check_setup.R` prints in its summary,
which reports READY, READY with a count of warnings, or NOT READY with the failing
stages listed.

### 17. Point of failure

Question: If you had trouble, where did it happen? Select all that apply.

Type: Choice, multiple answers, Other option on, optional

Subtitle: Leave blank if the setup completed without trouble.

Options: Installing R or Rtools / Installing the Xcode command line tools /
Installing R packages / Installing CmdStan / Running check_setup.R / No trouble

### 18. Setup check output

Question: Paste the summary block from the end of `check_setup.R` here.

Type: Text, long answer, optional

Subtitle: Everything from the Summary heading to the end. This identifies the
failing stage directly and is quicker to act on than a description of the error.
If the run reported NOT READY, email the whole output as well. Leave blank if you
have not run the check.

Added. `check_setup.R` already asks a participant to send its output back, and
the form is where that is received. Self-reported readiness in question 16 is
what a participant believes; this is what the machine reported.

## Handling the responses

Open the responses in Excel from the Responses tab, which writes a workbook to
OneDrive.

Responses are anonymous, so nobody can be contacted about their own answer. The
setup questions are read as counts across the room and answered by a reminder to
the whole list. Count the answers to question 16 that are not READY, and the
answers to question 15 that report installation needing IT approval. The pattern
in questions 17 and 18 says which stage is failing, which decides what that
reminder covers. An IT approval takes longer than anything else here, so send
that reminder first.

Where following up with an individual matters more than anonymity, add a short
answer question for an email address at the top of section 1 and mark it
optional. The Google Apps Script version in `scripts/build_participant_survey.gs`
still holds name and email as its first section.

Questions 1, 4, 9 and 10 together set how much time module 6 needs on priors and
how much module 2 needs before its first fit. Question 13 does the same for module
4.
