This is a python reconstruction of HM's initial Julia code (https://github.com/mberto79/LLT/blob/main/README.md).
Everything is formatted the same way, with the same variable names. 
The primary difference in python is the need to reference the functions created in function_definitions.py (using: fd."functionname"("variables")

example: 
# module imports
import function_definitions as fd  

# function call
C, D = fd.build_linear_system(y, span, Theta, c, m, a, a0)


