# AI Usage Report - Week 4

## Summary
Used Claude AI (Claude Haiku 4.5) and GPT 5/GPT 5 mini throughout Week 4 to implement and debug sequence alignment algorithms in both Python and Codon. Again as mentioned the last several weeks I find Claude does a much better job then other LLMs when dealing with Codon. The past couple weeks I ussed sonnet 4.5, but this week I used Anthropics newly released Haiku 4.5. I thought Sonnet 4.5 just realased as well but that's how fast this field works. 

Used Claude Haiku 4.5 for all the implementaitons here. 
## Algorithm implementation
**What I asked for**: Implementations of four sequence alignment algorithms: global, local, semi-global, and affine gap penalty alignment in Python.

**Claude provided**: 
Complete Python implementations. Was very good, with minimal errors.

## Codon conversion
**What I asked for**: (In the same claude chat) Convert x to codon ( where x = the above files)

**Claude provided**:
Complete Codon implementations. Again was very good with basically no errors. 

Used GPT 5/5 mini for the below as I like it better for things of explanatory nature. 
### Explaining
**What I asked for**: Help understanding why Python was running out of memory on large sequences(Pasted relevant info about the code and tests).

**What GPT explained**:
- O(n*m) space complexity means 16,569 × 16,499 = 273 million integers
- Python integer overhead (~28 bytes each) = ~7.6 GB per table
- Affine uses 3 tables = ~23 GB total, exceeding system RAM
