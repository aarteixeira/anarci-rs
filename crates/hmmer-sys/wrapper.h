/* bindgen entry point: pull in HMMER 3.4 + Easel public API.
 *
 * hmmer.h transitively includes the Easel headers it needs (esl_alphabet.h,
 * esl_sq.h, etc.) via p7_config.h's include path, but we include the Easel
 * headers we touch directly so their structs/functions (ESL_SQ, ESL_ALPHABET,
 * esl_sq_CreateFrom, esl_sq_Digitize, ...) are always in the bindgen output. */
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_getopts.h"
#include "hmmer.h"
