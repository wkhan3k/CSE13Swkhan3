# CSE13S-01 Winter 2026

## Assignment Directories

All work for an assignment must go into that assignment's corresponding folder. During
grading, only the contents of that folder will be checked, and anything outside of it
will not be graded. 

Not all of the assignment directories may be used (e.g. we may not have an assignment 8).
It's okay to have an empty directory for assignments that haven't been assigned yet. 

## Additional Files

There are a couple of extra files that we've provided for you.

`README.md`: That's the file you're reading now

`.gitignore`: This file tells git to ignore certain files, or certain files whose name
matches a particular pattern. For example, it's telling git to ignore all `.o` files,
meaning they won't show up when running `git status`, amongst other things. This helps
to ensure you don't accidentally commit files you shouldn't be committing. We've provided
some common extensions for git to ignore. If you wish, you can add more as you see fit.

`.clang-format`: This contains the formatting instructions used by `clang-format` so it 
knows how to format your code. Do not modify this file, as we require that your code
adheres to this format specification. 

You'll also find a blank README in each of the assignment directories. No need to worry
about those until you actually start each of those assignments. 

## Pipeline

Your repository is configured with an autograding system that you can run to check your 
work, known as the "pipeline". Instructions on how to use it will be provided during the
quarter. 

The pipeline makes use of certain secrets provided through Gitlab, which allows us to run
tests without revealing them. **DO NOT ATTEMPT TO MODIFY REPOSITORY OR PIPELINE SETTINGS 
TO GAIN ACCESS TO THESE SECRETS.** Any evidence of intentional tampering with the pipeline
to reveal these secrets will be treated as a
[Level 2 academic misconduct violation](https://ue.ucsc.edu/academic-misconduct.html), 
resulting in a 1--3 term suspension ***for a first offense***. 

There are multiple safeguards in place to prevent you from doing so, so you don't need
to worry about "accidentally" stumbling across the pipeline secrets. So long as you
use the pipeline according to directions, you will be fine. But woe unto you if you try
*anything* deliberate, which would be a clear and intentional violation of Sections
102.011, 102.014, 102.015, and 102.05 of the
[UCSC Policy on Student Conduct and Discipline](https://deanofstudents.ucsc.edu/student-conduct/student-handbook/100.003.pdf).

