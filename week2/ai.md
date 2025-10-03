# AI Usage Report

## AI Tools Used

### Claude Sonnet 4.5 (Primary)
Used extensively throughout the assignment for:

**Understanding the Codebase**
- Analyzing BioPython's motifs module structure
- Comparing working BioPython code with ported Codon versions
- Understanding MEME file format specifications
- Clarifying differences between counts and frequency matrices

**Debugging Type System Issues**
- Diagnosing Codon type inference problems with the Record class
- Troubleshooting `Dict[NoneType, NoneType]` vs `Dict[str, float]` 
type conflicts
- Attempting multiple solutions for Record class inheritance from 
list
- Understanding property initialization order and private 
attributes
- Debugging why `super().__init__()` vs `list.__init__(self)` made 
no difference

**Code Implementation Assistance**
- Porting background frequency parsing logic
- Implementing defensive parsing for optional MEME format 
parameters
- Converting Python duck typing to Codon's explicit type 
requirements
- Fixing alphabet handling (string vs list conversion)
- Correcting property setter usage vs direct method calls

**Git and Workflow Help**
- Removing embedded git repository (biopython folder)
- Creating .gitignore for compiled files
- Setting up evaluate.sh script
- Crafting CI workflow configuration
- Generating commit messages

**Documentation**
- Structuring this report with detailed technical gotchas
- Explaining each failed attempt to fix the Record class
- Documenting time estimates and lessons learned

### GPT-4 Mini (Minimal Usage)
Used briefly for:
- Quick syntax checks on markdown formatting
- Verifying bash script syntax for evaluate.sh

## What AI Could and Couldn't Do

**Successfully Helped With:**
- Identifying specific type system issues in error messages
- Suggesting multiple debugging approaches
- Explaining Codon vs Python differences
- Writing documentation and reports

**Could Not Solve:**
- The fundamental Record class type checking issue remains 
unsolved
- Multiple attempted fixes all failed despite AI suggestions
- This appears to be a genuine Codon type system limitation

## Estimate of AI vs Independent Work
- AI assistance: ~40% (debugging suggestions, documentation)
- Independent work: ~60% (actual code writing, testing, problem 
solving, trying solutions)

## Reflection
AI was helpful for understanding error messages and suggesting 
approaches, but the core type system issue proved unsolvable even 
with extensive AI assistance. This demonstrates that AI can 
accelerate debugging but cannot always overcome fundamental 
compiler limitations.
