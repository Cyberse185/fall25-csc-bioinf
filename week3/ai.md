# AI Usage Report

## Tools Used
- Claude (Sonnet 4.5)
- ChatGPT (GPT 5 mini)

As mentioned in the previous week, I find Sonnet 4.5 works much better for dealing with Codon than the few other LLM's I have tried. 

## How AI Was Used
- Understanding biotite's Cython source code structure
- Inital port of all the required code(sonnet 4.5: prompt: (assignment instructions) --- (the cython files))
- Debugging(sonnet 4.5 and gpt 5 mini) prompts were the error messages and the relevant sections of the code that the error messages mentioned. 
- For some reason, the models really struggled with what seemed like a easy error to fix. Recall from the report: Tried using custom exception classes but caused error: exceptions must derive from BaseException. For this error the models mostly just gave me garbage. 

