#Program Name: Strassen.py
#Program Desc: To Find Product of Matrices using Strassen Method
#Author: Sunny Kumar
#Reg=201700282
#Section: B
#Roll NO: 13
#
#Implemented for size=2^n(i.e:1,2,4,8,16,32,...)
#
import time     #To Measure the Time
import random   #To Create Random numbers as the input
import sys      #To increase the Recursion Limit


sys.setrecursionlimit(1000000)  #increasing the Recursion limit


def MatMul(A,B,n):  #Trivial Matrix Multiplication Method
    C=[]
    for i in range(n):
        temp=[]
        for j in range(n):
            sm=0
            for k in range(n):
                sm+=A[i][k]*B[k][j]
            temp.append(sm)
        C.append(temp)
    return C

def Allocate(n):
    return[[0 for i in range(n)] for j in range(n)]

#For declaring the Quardrants
def Quad(A,B,size): 
    n=size//2
    A11 = Allocate(n)
    A12 = Allocate(n)
    A21 = Allocate(n)
    A22 = Allocate(n)
    B11 = Allocate(n)
    B12 = Allocate(n)
    B21 = Allocate(n)
    B22 = Allocate(n)

    #Assigning the Values of the Quardrants
    for i in range(n):
        for j in range( n):
            A11[i][j] = A[i][j]         
            A12[i][j] = A[i][j + n]     
            A21[i][j] = A[i + n][j]     
            A22[i][j] = A[i + n][j + n] 

            B11[i][j] = B[i][j]          
            B12[i][j] = B[i][j + n]    
            B21[i][j] = B[i + n][j]   
            B22[i][j] = B[i + n][j + n] 
    
    return(A11,A12,A21,A22,B11,B12,B21,B22)

#Matrix Addition
def addMat(A,B):    
    n=len(A)
    C = Allocate(n)
    for i in range(n):
        for j in range(n):
            C[i][j]=A[i][j]+B[i][j]
    return C

#Matrix Subtraction
def subMat(A,B):
    n=len(A)
    C = Allocate(n)
    for i in range(n):
        for j in range(n):
            C[i][j]=A[i][j]-B[i][j]
    return C

#Combing the Matrices C11...C22
def Combine(C11,C12,C21,C22,n):
    C = Allocate(n)
    k=int(n/2)
    for i in range(k):
        for j in range(k):
            C[i][j] = C11[i][j]
            C[i][j + k] = C12[i][j]
            C[i + k][j] = C21[i][j]
            C[i + k][j + k] = C22[i][j]
    return C

#Func Strassen will divide the Matrices 7 times
def Strassen(A,B):
    n=len(A)
    C=Allocate(n)
    if n==1:                    #only one element is in the Matrices A and B
        C[0][0]=A[0][0]*B[0][0] 
    else: 
        #Dividing the Matrices Indices by 2
        mid=int(n/2)

        #Creating The Smaller Quadrants of the Matrices
        A11,A12,A21,A22,B11,B12,B21,B22=Quad(A,B,n)
        #Creating SubSection of the Product Matrix
        C11 = Allocate(mid)
        C12 = Allocate(mid)
        C21 = Allocate(mid)
        C22 = Allocate(mid)

        #Creating matrices S1,S2......S10
        S1 = Allocate(mid)
        S2 = Allocate(mid)
        S3 = Allocate(mid)
        S4 = Allocate(mid)
        S5 = Allocate(mid)
        S6 = Allocate(mid)
        S7 = Allocate(mid)
        S8 = Allocate(mid)
        S9 = Allocate(mid)
        S10 = Allocate(mid)

        #Allocating Matrices P1...P7 for Strassen Multiplication
        P1 = Allocate(mid)
        P2 = Allocate(mid)
        P3 = Allocate(mid)
        P4 = Allocate(mid)
        P5 = Allocate(mid)
        P6 = Allocate(mid)
        P7 = Allocate(mid)
        
        #Assigning the Values of S1...S10 
        S1=subMat(B12,B22)
        S2=addMat(A11,A12)
        S3=addMat(A21,A22)
        S4=subMat(B21,B11)
        S5=addMat(A11,A22)
        S6=addMat(B11,B22)
        S7=subMat(A12,A22)
        S8=addMat(B21,B22)
        S9=subMat(A11,A21)
        S10=addMat(B11,B12)
        
        #Recursive Calls for breaking the Matrices into smaller groups
        P1=Strassen(A11,S1)
        P2=Strassen(S2,B22)
        P3=Strassen(S3,B11)
        P4=Strassen(A22,S4)
        P5=Strassen(S5,S6)
        P6=Strassen(S7,S8)
        P7=Strassen(S9,S10)

        #Assinging the Values of the SubMatrices C11...C22
        C12=addMat(P1,P2)
        C11=subMat((addMat(addMat(P4,P5),P6)), P2)
        C21=addMat(P3,P4)
        C22=subMat((addMat(P1,P5)),addMat(P3,P7))

        #Merging the Matrices C11...C22 to C
        C=Combine(C11,C12,C21,C22,n)
    return C


#Start of the Programme    
print("**********\n\nImplemented for size=2^n,(2,4,8,...)\n\n**********\n")
tC=int(input("Enter the Number of Test Cases:"))


StMM=[] #Stores the Execution Time for Strassen Method
MM=[]   #Stores the Execution Time for Trivial Matrix Multiplication

for t in range(tC):
    print("Test Case ",t+1)
    n=int(input("Enter the Dimension of the Square Matrix: "))
    A = Allocate(n)
    for i in range(n):
        for j in range(n):
            A[i][j]=random.randrange(100)
  
    B = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            B[i][j]=random.randrange(100)
    #Calculation of time for Strassen Method   
    st=time.time()
    M1=Strassen(A,B)
    ed=time.time()
    StMM.append(ed-st)
    #Calculation of time for Trivial Method
    st=time.time()
    M2=MatMul(A,B,n)
    ed=time.time()
    MM.append(ed-st)
print("\n\nTime Consumed Strassen Method: ",StMM)
print("\n\nTime Consumed Trivial Method: ",MM)
