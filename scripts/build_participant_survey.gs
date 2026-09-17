/**
 * Builds the pre-workshop participant survey as a Google Form.
 *
 * Course: concentration-response modelling, SETAC Australasia, September 2026.
 * Presenter: Dr Rebecca Fisher, Australian Institute of Marine Science.
 *
 * To run: open script.google.com, create a new project, replace the contents of
 * Code.gs with this file, then run buildSurvey(). Authorise when prompted. The
 * form and its response spreadsheet are created in the Drive root of whichever
 * account ran it, and the URLs are printed to the execution log (View >
 * Executions, or Ctrl+Enter).
 *
 * Running it twice creates a second form. To change wording after the fact,
 * edit the form in the Forms editor rather than re-running this.
 */

var CONFIG = {
  title: 'Concentration-response modelling workshop: participant survey',

  // Items 9 to 13 of the question set: short checks on concepts the course
  // teaches. Set to false to issue the shorter background-only form and run
  // these live on the day instead.
  includeConceptChecks: true,

  // Creates a spreadsheet and links responses to it. Without this the responses
  // live inside the form and are exported by hand.
  createResponseSheet: true,

  // Shown after submission. No congratulation, no exclamation.
  confirmation: 'Your response has been recorded. If your setup check did not ' +
    'report READY, we will contact you before the workshop.'
};

function buildSurvey() {
  var form = FormApp.create(CONFIG.title);

  form.setDescription(
    'This survey shapes the pacing and emphasis of the workshop. It takes about ' +
    'four minutes. The software setup section at the end matters most: a ' +
    'participant whose C++ toolchain is not working cannot fit models on the ' +
    'day, and there is no troubleshooting slot in the programme.');

  form.setProgressBar(true);
  form.setAllowResponseEdits(true);
  form.setConfirmationMessage(CONFIG.confirmation);

  // --- Identification --------------------------------------------------------
  // An explicit email question rather than form.setCollectEmail(), which
  // requires a Google sign-in and excludes anyone without an account.
  form.addTextItem()
    .setTitle('Name')
    .setRequired(true);

  form.addTextItem()
    .setTitle('Email address')
    .setHelpText('Used to contact you before the workshop if your software ' +
                 'setup needs attention.')
    .setValidation(FormApp.createTextValidation()
      .setHelpText('Enter a valid email address.')
      .requireTextIsEmail()
      .build())
    .setRequired(true);

  // --- Background and prior experience --------------------------------------
  form.addPageBreakItem()
    .setTitle('Background and prior experience');

  mc(form, 'Which best describes your current experience with R?', [
    'Never used it',
    'Use it occasionally',
    'Use it regularly',
    'Expert'
  ], { required: true });

  mc(form, 'Have you fitted a concentration-response or dose-response model ' +
           'before, in any software?', [
    'Never',
    'A few times',
    'Routinely'
  ], { required: true });

  cb(form, 'Which of these have you used before? Select all that apply.', [
    'drc',
    'bayesnec',
    'brms',
    'Stan directly, through rstan or cmdstanr',
    'None of these'
  ], { other: true, required: true,
       help: 'Use Other to name software outside R, for example ToxRat, CETIS ' +
             'or GraphPad.' });

  mc(form, 'How would you describe your comfort with Bayesian statistics?', [
    'New to me',
    'Know the basic idea',
    'Have interpreted posteriors before',
    'Fit Bayesian models routinely'
  ], { required: true });

  mc(form, 'What is your day-to-day role?', [
    'Regulatory risk assessor',
    'Consultant',
    'Academic researcher',
    'Lab-based ecotoxicologist'
  ], { other: true, required: true });

  // --- Goals for the day -----------------------------------------------------
  form.addPageBreakItem()
    .setTitle('Goals for the day');

  cb(form, 'What are you most hoping to get from the workshop?', [
    'Confidence installing and running bayesnec myself',
    'Understanding when NEC and when ECx is appropriate',
    'Justifying a model choice to a reviewer or regulator',
    'Comparing Bayesian and frequentist results',
    'Applying the methods to my own data'
  ], { other: true, required: true, atMost: 2,
       help: 'Select up to two.' });

  mc(form, 'Do you currently need to report NOEC, ECx, or NEC and NSEC values ' +
           'in your work?', [
    'NOEC only',
    'ECx',
    'NEC or NSEC',
    'A mix, depending on the study',
    'Not sure what these are'
  ], { required: true });

  form.addParagraphTextItem()
    .setTitle('Is there a specific dataset or problem you are hoping to apply ' +
              'this to afterwards?')
    .setHelpText('A sentence is enough. Response type, number of ' +
                 'concentrations, and what you need to estimate.');

  // --- Concept checks --------------------------------------------------------
  if (CONFIG.includeConceptChecks) {
    form.addPageBreakItem()
      .setTitle('Concept checks')
      .setHelpText('These set the starting point for the relevant modules. ' +
                   'Answer as you stand now, before the course covers them. ' +
                   'All optional.');

    mc(form, 'Do you know what distinguishes a NEC model from an ECx model?', [
      'Yes',
      'Roughly',
      'No'
    ]);

    mc(form, 'How comfortable are you reading a posterior distribution or a ' +
             'credible interval?', [
      'Very',
      'Somewhat',
      'Not really'
    ]);

    mc(form, 'Have you had to justify a choice of statistical distribution, ' +
             'such as binomial, beta or Gaussian, for response data?', [
      'Yes, routinely',
      'Once or twice',
      'Never'
    ]);

    mc(form, 'Do you have a view on informative against weakly informative ' +
             'priors?', [
      'Have a view',
      'Aware of the question',
      'New territory'
    ]);

    mc(form, 'How important is model averaging, against selecting a single ' +
             'best model, in your own work?', [
      'Essential',
      'Useful',
      'Not something I currently need',
      'Not sure what model averaging is'
    ]);
  }

  // --- Software setup --------------------------------------------------------
  form.addPageBreakItem()
    .setTitle('Software setup')
    .setHelpText('The setup instructions and check_setup.R were issued with ' +
                 'the joining email. Run the check on the machine you will ' +
                 'bring on the day.');

  mc(form, 'Which operating system is on the machine you will bring?', [
    'Windows',
    'macOS, Apple silicon',
    'macOS, Intel',
    'Linux'
  ], { other: true, required: true,
       help: 'The toolchain differs by platform: Rtools on Windows, the Xcode ' +
             'command line tools on macOS.' });

  mc(form, 'Can you install software on that machine without an IT request?', [
    'Yes',
    'No, installation needs IT approval',
    'Not sure'
  ], { required: true,
       help: 'A managed laptop that refuses a compiler install is the most ' +
             'common reason a participant cannot fit a model on the day. ' +
             'Raise the request now if the answer is no.' });

  mc(form, 'Did check_setup.R report READY?', [
    'Yes, READY',
    'READY, with warnings',
    'NOT READY',
    'Attempted, but the script would not run',
    'Have not tried yet'
  ], { required: true });

  cb(form, 'If you had trouble, where did it happen? Select all that apply.', [
    'Installing R or Rtools',
    'Installing the Xcode command line tools',
    'Installing R packages',
    'Installing CmdStan',
    'Running check_setup.R',
    'No trouble'
  ], { other: true,
       help: 'Leave blank if the setup completed without trouble.' });

  form.addParagraphTextItem()
    .setTitle('Paste the output of check_setup.R here.')
    .setHelpText('The whole output, from the R version line to the final ' +
                 'summary. This is the single most useful thing you can give ' +
                 'us: it identifies the failing stage directly and is far ' +
                 'quicker to act on than a description of the error. Leave ' +
                 'blank if you have not run it.');

  // --- Response destination --------------------------------------------------
  if (CONFIG.createResponseSheet) {
    var sheet = SpreadsheetApp.create(CONFIG.title + ' (responses)');
    form.setDestination(FormApp.DestinationType.SPREADSHEET, sheet.getId());
    Logger.log('Responses sheet : ' + sheet.getUrl());
  }

  Logger.log('Edit the form    : ' + form.getEditUrl());
  Logger.log('Send this link   : ' + form.shortenFormUrl(form.getPublishedUrl()));
}

/**
 * Adds a single-answer question. opts: required, other, help.
 */
function mc(form, title, choices, opts) {
  opts = opts || {};
  var item = form.addMultipleChoiceItem()
    .setTitle(title)
    .setChoiceValues(choices)
    .setRequired(opts.required === true);
  if (opts.other) item.showOtherOption(true);
  if (opts.help) item.setHelpText(opts.help);
  return item;
}

/**
 * Adds a multiple-answer question. opts: required, other, help, atMost.
 */
function cb(form, title, choices, opts) {
  opts = opts || {};
  var item = form.addCheckboxItem()
    .setTitle(title)
    .setChoiceValues(choices)
    .setRequired(opts.required === true);
  if (opts.other) item.showOtherOption(true);
  if (opts.help) item.setHelpText(opts.help);
  if (opts.atMost) {
    item.setValidation(FormApp.createCheckboxValidation()
      .setHelpText('Select at most ' + opts.atMost + '.')
      .requireSelectAtMost(opts.atMost)
      .build());
  }
  return item;
}
