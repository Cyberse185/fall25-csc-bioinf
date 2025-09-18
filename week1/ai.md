# AI Usage Documentation

## AI Tools Used
- **Claude (Anthropic)** – Helped with Codon implementation, debugging, and workflow issues  
  - Version: Claude 3.5 Sonnet (accessed via claude.ai)
- **ChatGPT (OpenAI)** – Helped with Codon implementation, `evaluate.sh`, and CI 
troubleshooting  
  - Version: GPT-5 min (accessed via chat.openai.com)

## Key Areas of AI Assistance

### Codon Implementation
- Asked about translating Python code to Codon  
- Pasted code snippets that weren’t working and got guidance to fix syntax and imports  
- Debugged runtime issues and library setup for Codon  

### Automation Script (`evaluate.sh`)
- Asked how to run Python and Codon automatically and compare outputs  
- Asked how to format results as a table for CI  

### Debugging and CI/CD
- Asked why GitHub Actions failed to find files or dependencies  
- Asked how to install Python dependencies and set CODON_PYTHON on the CI runner  

### Repository Organization and Documentation
- Asked whether files were in the correct structure for CI  

## Sample Prompts Used
- "I'm getting a libpython error when running Codon, how do I fix this?"
- "My Codon code compiles but segfaults on large data, what could be the issue?"
- "Help me implement `evaluate.sh` to run Python and Codon automatically and print N50 table" — 
GPT  
- "Why is GitHub Actions failing on `cd week1` or missing dependencies?" — Claude  
- "How do I fix Codon import and syntax errors?(pasted code)" — GPT  
