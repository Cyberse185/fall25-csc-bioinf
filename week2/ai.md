# AI Usage Documentation

## AI Tools Used

### Claude Sonnet 4.5 (Primary Tool)
* **Version**: Claude Sonnet 4.5 (newly released - accessed via claude.ai)
* **Primary Usage**: 80-85% of AI assistance
* **Performance**: Found this model significantly more effective than previous versions for understanding Codon's type system and debugging complex compilation errors

### ChatGPT (Secondary Tool)
* **Version**: GPT-4o mini (accessed via chat.openai.com)
* **Secondary Usage**: 15-20% of AI assistance
* **Use Cases**: Second opinion

## Claude Sonnet 4.5 Performed Better

Main thing to take away here is that I found Claude Sonnet 4.5 was far superior than other models for this task. 

## Key Areas of AI Assistance

### Initial Code Generation
* Asked Claude to generate skeleton files for all BioPython ports
* Requested basic class structures for Motif, Record, and matrix classes
* Generated initial property getter/setter templates
* Created starting point for minimal.py parser structure

**Sample Prompts**:
* "Generate a skeleton Motif class for Codon based on BioPython's structure"
* "Create basic FrequencyPositionMatrix class with normalize and log_odds methods"
* "Help me set up the initial structure for minimal.py parser"

### Understanding BioPython Codebase
* Asked Claude to analyze BioPython's motifs module structure
* Compared working BioPython __init__.py with my ported version
* Asked about differences between frequency matrices and count matrices
* Requested explanation of MEME file format specifications

**Sample Prompts**:
* "I attached the working BioPython minimal.py and my non-working version, what are the key differences?"
* "Why does BioPython's FrequencyPositionMatrix store counts instead of frequencies?"

### Debugging Codon Type System Issues
* Pasted full error messages and asked for interpretation
* Shared my Record class code that failed type checking
* Asked why `Dict[NoneType, NoneType]` was being inferred
* Requested explanations of property initialization order
* Iteratively debugged through 6+ different attempted fixes

**Sample Prompts**:
* "Getting 'cannot typecheck class Record(list):' error - here's my code (pasted)"
* "Tried adding type annotations but still getting Dict[NoneType, NoneType], what's wrong?"
* "Why does `self.background = {}` infer the wrong type even with explicit annotations?"
* "That didn't work, getting same error - what else can I try?"

### Code Implementation Help
* Asked how to port background frequency parsing from BioPython
* Requested help converting setter methods to property decorators
* Asked about defensive parsing for optional MEME parameters
* Debugged alphabet handling (string vs list conversion)
* Fixed property vs method call issues

**Sample Prompts**:
* "How do I parse 'Background letter frequencies' followed by multiple lines of data?"
* "Should I use `motif.background =` or `motif.set_background()`?"
* "Need to handle optional 'nsites=' parameter that might not exist in line"
* "Getting AttributeError on _pseudocounts - help debug initialization order"

### Git and CI Workflow
* Asked how to remove embedded git repository (biopython folder)
* Requested help creating .gitignore for compiled files
* Asked about setting up evaluate.sh script
* Debugged GitHub Actions workflow configuration

**Sample Prompts** (mostly to GPT-4o mini):
* "Getting 'adding embedded git repository' warning, how to fix?"
* "Need bash script to run both Python and Codon tests"
* "How to add BioPython to CI dependencies?"

## What AI Could and Couldn't Do

### Successfully Helped With:
* Generating initial skeleton code for all ported modules
* Interpreting cryptic Codon error messages
* Suggesting multiple debugging approaches (even though most failed)
* Explaining fundamental differences between Python and Codon type systems
* Writing clear documentation of technical issues
* Git workflow troubleshooting
* Iterative debugging through many failed attempts

### Could Not Solve:
* The Record class type checking issue remains unsolved after 6+ different approaches
