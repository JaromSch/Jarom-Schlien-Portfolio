from sympy import symbols, cos, sin, pi
import itertools
import math
m, n, p = symbols('m n p', real=True, nonnegative=True)


#---- Assumption: Symmetric and Balanced Stack ----#


#---- Type in: What number of plies do you want to investigate? ----#
max_plies = 8


#---- These are example values for material, ply thickness, width and length ----#
# Declaring variables; R corresponds to beta (stress ratio) and u corresponds to alpha (aspect ratio)
# m & p = number of half waves in x and y, n = number of plies
a = 400
b = 150
t = 0.2
Nx = -0.5
Ny = 0.2
Tau = 0.1
R = Ny / Nx
Q11 = 131098
Q12 = 3327
Q22 = 10084
Q66 = 5000
u = a / b


# Since we have a fixed material, we use general variables for the Q matrix. I didn't implement the solver for the Q-matrix here. The analyst would therefore have to do that manually

#Q-bar for zero degrees. Fixed values dependent on material and orientation
Q110 = Q11 * (cos(0) ** 4) + 2 * (Q12 + 2 * Q66) * (sin(0) ** 2) * (cos(0) ** 2) + Q22 * (sin(0) ** 4)
Q120 = (Q11 + Q22 - 4 * Q66) * (sin(0) ** 2) * (cos(0) ** 2) + Q12 * ((sin(0) ** 4) + (cos(0) ** 4))
Q220 = Q11 * (sin(0) ** 4) + 2 * (Q12 + 2 * Q66) * (sin(0) ** 2) * (cos(0) ** 2) + Q22 * (cos(0) ** 4)
Q660= (Q11 + Q22 - 2 * Q12 - 2 * Q66) * (sin(0) ** 2) * (cos(0) ** 2) + Q66 * ((sin(0) ** 4) + cos(0) ** 4)

#Q-bar for 45 degrees. Fixed values dependent on material and orientatio
Q1145 = Q11 * (cos(pi/4) ** 4) + 2 * (Q12 + 2 * Q66) * (sin(pi/4) ** 2) * (cos(pi/4) ** 2) + Q22 * (sin(pi/4) ** 4)
Q1245 = (Q11 + Q22 - 4 * Q66) * (sin(pi/4) ** 2) * (cos(pi/4) ** 2) + Q12 * ((sin(pi/4) ** 4) + (cos(pi/4) ** 4))
Q2245 = Q11 * (sin(pi/4) ** 4) + 2 * (Q12 + 2 * Q66) * (sin(pi/4) ** 2) * (cos(pi/4) ** 2) + Q22 * (cos(pi/4) ** 4)
Q6645= (Q11 + Q22 - 2 * Q12 - 2 * Q66) * (sin(pi/4) ** 2) * (cos(pi/4) ** 2) + Q66 * ((sin(pi/4) ** 4) + cos(pi/4) ** 4)

#Q-bar for -45 degrees. Fixed values dependent on material and orientation
Q1145n = Q11 * (cos(-pi/4) ** 4) + 2 * (Q12 + 2 * Q66) * (sin(-pi/4) ** 2) * (cos(-pi/4) ** 2) + Q22 * (sin(-pi/4) ** 4)
Q1245n = (Q11 + Q22 - 4 * Q66) * (sin(-pi/4) ** 2) * (cos(-pi/4) ** 2) + Q12 * ((sin(-pi/4) ** 4) + (cos(-pi/4) ** 4))
Q2245n = Q11 * (sin(-pi/4) ** 4) + 2 * (Q12 + 2 * Q66) * (sin(-pi/4) ** 2) * (cos(-pi/4) ** 2) + Q22 * (cos(-pi/4) ** 4)
Q6645n= (Q11 + Q22 - 2 * Q12 - 2 * Q66) * (sin(-pi/4) ** 2) * (cos(-pi/4) ** 2) + Q66 * ((sin(-pi/4) ** 4) + cos(-pi/4) ** 4)

#Q-bar for 90 degrees. Fixed values dependent on material and orientation
Q1190 = Q11 * (cos(pi/2) ** 4) + 2 * (Q12 + 2 * Q66) * (sin(pi/2) ** 2) * (cos(pi/2) ** 2) + Q22 * (sin(pi/2) ** 4)
Q1290 = (Q11 + Q22 - 4 * Q66) * (sin(pi/2) ** 2) * (cos(pi/2) ** 2) + Q12 * ((sin(pi/2) ** 4) + (cos(pi/2) ** 4))
Q2290 = Q11 * (sin(pi/2) ** 4) + 2 * (Q12 + 2 * Q66) * (sin(pi/2) ** 2) * (cos(pi/2) ** 2) + Q22 * (cos(pi/2) ** 4)
Q6690= (Q11 + Q22 - 2 * Q12 - 2 * Q66) * (sin(pi/2) ** 2) * (cos(pi/2) ** 2) + Q66 * ((sin(pi/2) ** 4) + cos(pi/2) ** 4)

#it probably would have been more efficient to initialize four calculations and iterate with the respective angle to get all the Q values. Next time :)

#explanation: Within the code I added lines that are meant to explain certain lines.
#Idea: Due to practical reasons, I only analyze plies with the orientations 0°, +-45° and 90°. However, adding other orientations wouldn't be a big additional coding effort.

#In the first step (def generate_balanced_symmetric_stacks), all possible stacking sequences are computed and put into a list.
    #According to the constraints balanced, symmetric and minimum 10% ply share, the list gets filtered in the same step
#In the second step, all relevant D-values get computed which are D11, D12, D22 and D66. I named them just D1, D2, D3 and D4 for simplicity. Same notation is used for Q-values.
    #Into the computation I added a very important step: As I don't want to analyze the critical stress but rather the stress/weight ratio, I divide every D-value by the number of plies.
    #Weight = density * volume = height * width * length * density. We assume constant width, length, thickness and unit density. Height = thickness * number of plies
    #Mathematically, dividing every D-value by n is the same as dividing the final stress by n
    #Depending on the position of the respective ply, the correct Q-value gets chosen from the paragraph/code below
#In the third step we perform the actual stress calculations
    #Second paragraph: formular for critical stress ("compute_stress")
    #First paragraph: Take the result and iterate to find m and p over a 4x4 table to get the minimal maximum stress, which is what we want at the end
    #Third paragraph: Give "max_stress" the value 0. Now we iterate over all combinations and overwrite "max_stress" to achieve the highest critical stress


#1st step: Computing possible sequences and filtering
def generate_balanced_symmetric_stacks(max_plies):
    orientations = [0, 45, -45, 90]
    all_stacks = []

    for n in range(1, max_plies + 1):
        for comb in itertools.product(orientations, repeat=n):
            stack = list(comb)
    #You can add print(stack) here to see all possible combinations
            if is_balanced(stack) and is_symmetric(stack) and satisfies_minimum_occurrence(stack, n):
    #You can add print(stack) here to see all valid combinations
                all_stacks.append(stack)
    return all_stacks
    #balanced condition
def is_balanced(stack):
    return stack.count(45) == stack.count(-45)
    #symmetric condition
def is_symmetric(stack):
    return stack == stack[::-1]
    #10 percent ply share condition
def satisfies_minimum_occurrence(stack, n):
    min_count = max(1, n // 10)
    for orientation in set(stack):
        if stack.count(orientation) < min_count:
            return False
    return True


#2nd step: Calculating D-values
def calculate_D_values(stack, t):
    n = len(stack)
    D1, D2, D3, D4 = 0, 0, 0, 0
    # for readability: D11=D1 D12=D2 D22=D3 D66=D4

    for k in range(1, n + 1):
    #k ranges from 1 to n with n being the number of plies
        Q_values = get_Q_values(stack[k - 1])
    #Get the Q values for the respective ply at position k
        term = ((k * t - n * t * 0.5) ** 3 - ((k - 1) * t - n * t * 0.5) ** 3) / (3*n)
    #The magic happens when dividing by 3*n instead of just 3 at the very end of "term" to get stress/weight ratio
    # Calculating each D value
        D1 += Q_values['Q1'] * term
        D2 += Q_values['Q2'] * term
        D3 += Q_values['Q3'] * term
        D4 += Q_values['Q4'] * term
    return D1, D2, D3, D4


    #That's where "calculate_D_values" gets the respective Q-values from
def get_Q_values(orientation):

    if orientation == 0:
        return {'Q1': Q110, 'Q2': Q120, 'Q3': Q220, 'Q4': Q660}
    elif orientation == 45:
        return {'Q1': Q1145, 'Q2': Q1245, 'Q3': Q2245, 'Q4': Q6645}
    elif orientation == -45:
        return {'Q1': Q1145n, 'Q2': Q1245n, 'Q3': Q2245n, 'Q4': Q6645n}
    elif orientation == 90:
        return {'Q1': Q1190, 'Q2': Q1290, 'Q3': Q2290, 'Q4': Q6690}


#3rd step: perform critical stress calculation

    #minimal maximum critical stress
def calculate_critical_stress(D1, D2, D3, D4, m, p):
    min_stress = float('inf')
    #defining the minimal stress as an arbitrary high number

    #Trial and error to get minimum value of half waves m and p. Iterating through a 4x4 table
    for m in range(1, 5):
        for p in range(1, 5):

            stress = compute_stress(D1, D2, D3, D4, m, p)

            #Overwrite the minimum stress if new computed stress is lower => finding correct m and p
            if stress < min_stress:
                min_stress = stress
    return min_stress
import math
    #This is the actual formular for calculating the critical stress
def compute_stress(D1, D2, D3, D4, m, p):

    critical_stress = (math.pi**2)*((D1*(m/u)**4)+2*(D2+D4)*((m*p/u)**2)+D3*p**4)/(((t*(b**2))*(((m/u)**2)+R*(p**2))))

    return abs(critical_stress)

    #Here, the final critical stress is found by comparing all iterations and overwriting smaller values
def find_max_stress(max_plies, t):
    max_stress = 0
    best_stack = None

    all_stacks = generate_balanced_symmetric_stacks(max_plies)

    for stack in all_stacks:
        D1, D2, D3, D4 = calculate_D_values(stack, t)
        stress = calculate_critical_stress(D1, D2, D3, D4, m, p)

        if stress > max_stress:
            max_stress = stress
            best_stack = stack

    return max_stress, best_stack


max_stress, best_stack = find_max_stress(max_plies, t)
print()
print("Maximum number of plies checked:", max_plies)
print("Number of plies leading to the best stress-mass ratio:", len(best_stack))
print()
print("Critical stress-mass ratio:", max_stress)
print("Best stacking sequence:", best_stack)
print()
print("Have a safe flight!")



