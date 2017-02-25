/*
	A tool that contains numerical Methods
	@author
	Name			: pleycpl
	go version 		: go1.5.4 linux/amd64
	gcc --version   : gcc (GCC) 5.3.1 20160406 (Red Hat 5.3.1-6)
*/
package main

import (
	"fmt"
	"os"
)

type Necessary struct { //General struct for Methods
	equationSize     int
	equationElements [10]int
	x0               float32
	x0new            float32
	x1               float32
	x2               float32
	deltaX           float32
	epsilon          float32
	firstRange       float32
	secondRange      float32
	thirdRange       float32
}
type Matrix struct { // Matrix struct for Gaussian(Jordan)Elimination Methods
	matrixElements [10][10]float32
	matrixSize     int
	arrayElements  [10]float32
}

//Methods
func (itDoesNotMatter Matrix) determinat() (det float32) {
	A := itDoesNotMatter.matrixElements
	size := itDoesNotMatter.matrixSize
	itDoesNotMatter.displayMatrix()
	if size == 1 {
		det = A[0][0]
	} else if size == 2 {
		det = A[1][1]*A[0][0] - A[1][0]*A[0][1]
	} else { //matrix 3x3, 4x4, 5x5 ....
		B := Matrix{}
		B.displayMatrix()
		C := 1
		for i := 0; i < size; i++ {
			for j := 1; j < size; j++ {
				for k := 0; k < size-1; {
					if k != i {
						B.matrixElements[j-1][k] = A[j][i]
						k++
					}
					B.displayMatrix()
				}
			}
			B.matrixSize = size - 1
			det += float32(C) * A[i][i] * B.determinat() //Recursive Methods
			C *= -1
		}
	}
	return det
	//Link:https://www.goinggo.net/2013/09/recursion-and-tail-calls-in-go_26.html
}
func (itDoesNotMatter Matrix) displayArray() {
	fmt.Println("*** Array Elements ***")
	for i := 0; i < itDoesNotMatter.matrixSize; i++ {
		fmt.Println(itDoesNotMatter.arrayElements[i])
	}
}
func (itDoesNotMatter Matrix) displayMatrix() {
	fmt.Println("*** Matrix Elements ***")
	for i := 0; i < itDoesNotMatter.matrixSize; i++ {
		for j := 0; j < itDoesNotMatter.matrixSize; j++ {
			fmt.Printf("%3f ", itDoesNotMatter.matrixElements[i][j])
		}
		fmt.Println()
	}
}
func (itDoesNotMatter *Matrix) takeArrayElements() {
	fmt.Println("Please Entry to Elements of Array")
	for i := 0; i < itDoesNotMatter.matrixSize; i++ {
		fmt.Printf("Array[%d] : ", i)
		fmt.Scanln(&itDoesNotMatter.arrayElements[i])
	}
}
func (itDoesNotMatter *Matrix) takeMatrixElements() {
	fmt.Println("Please Entry to Elements of Matrix...")
	for i := 0; i < itDoesNotMatter.matrixSize; i++ {
		for j := 0; j < itDoesNotMatter.matrixSize; j++ {
			fmt.Printf("Matrix[%d][%d] : ", i, j)
			fmt.Scanln(&itDoesNotMatter.matrixElements[i][j])
		}
	}
}
func (itDoesNotMatter *Matrix) takeMatrixSize() { //Denklemleriniz kac tane bilinmeyenden olusuyor.
	fmt.Println("Please Entry to Size of Matrix: ")
	fmt.Scan(&itDoesNotMatter.matrixSize)
}
func (itDoesNotMatter *Necessary) takeRanges() {
	fmt.Print("Please Entry to First Range : ")
	fmt.Scan(&itDoesNotMatter.firstRange)
	fmt.Print("Please Entry to Second Range : ")
	fmt.Scan(&itDoesNotMatter.secondRange)
}
func (itDoesNotMatter *Necessary) takeEpsilon() {
	fmt.Print("Please Entry to Epsilon : ")
	fmt.Scan(&itDoesNotMatter.epsilon)
}
func (itDoesNotMatter *Necessary) takeDeltaX() {
	fmt.Print("Please Entry to DeltaX : ")
	fmt.Scan(&itDoesNotMatter.deltaX)
}

/*
func (itDoesNotMatter *Necessary) takeX2() {
	fmt.Print("Please Entry to x2 : ")
	fmt.Scan(&itDoesNotMatter.x2)
}
*/
func (itDoesNotMatter *Necessary) takeX1() {
	fmt.Print("Please Entry to x1 : ")
	fmt.Scan(&itDoesNotMatter.x1)
}
func (itDoesNotMatter *Necessary) takeX0() {
	fmt.Print("Please Entry to x0 : ")
	fmt.Scan(&itDoesNotMatter.x0)
}
func (itDoesNotMatter *Necessary) takeSize() {
	fmt.Print("Please Entry to Equation Size : ")
	fmt.Scan(&itDoesNotMatter.equationSize)
}
func (itDoesNotMatter *Necessary) takeElements() {
	fmt.Println("Please Entry to Equation Elements,but left to right...")
	for i := 0; i < itDoesNotMatter.equationSize; i++ {
		fmt.Scan(&itDoesNotMatter.equationElements[i])
	}
}

//Functions
func main() {
	fmt.Print("\033[2J") //ClearTheScreenWith ANSI escape code Source:http://rosettacode.org/wiki/Terminal_control/Clear_the_screen#Go
	switch tmp := Menu(); tmp {
	case 1:
		GraphicalAnalysis()
	case 2:
		BiSection()
	case 3:
		RegulaFalsi()
	case 4:
		NewtonRaphson()
	case 5:
		Secant()
	case 6:
		GaussianJordanElimination()
	case 7:
		GaussianElimination()
	case 9:
		Trapez()
	case 10:
		Simpson()
	case 11:
		CentralDifference()
	case 0:
		os.Exit(0)
	default:
		fmt.Println("Please Entry between 1-10 ")
	}
}

// Supporting Functions
func fxCalculate(equation Necessary, value float32) float32 { //Calculate to Fx equation from value
	var result float32
	var power float32 = 1
	for i := 0; i < equation.equationSize; i++ {
		result += float32(equation.equationElements[i]) * power
		power *= value
	}
	return result
}
func Abs(x float32) float32 { //Absolute value . if - then + else +
	if x < 0 {
		return -x
	}
	if x == 0 {
		return 0
	}
	return x
}
func derivativeCalculate(equation Necessary, value float32) float32 { // Calculate to Fx equation f
	return (fxCalculate(equation, value) - fxCalculate(equation, value-0.001)) / 0.001
}

// Main Functions

func CentralDifference() { //Numerical Derivative Method
	p := Necessary{}
	p.takeSize()
	p.takeElements()
	p.takeX0()
	p.takeEpsilon() // difference is epsilon
	solutionWithMethods := (fxCalculate(p, p.x0+p.epsilon/2) - fxCalculate(p, p.x0-p.epsilon/2)) / p.epsilon
	realSolution := derivativeCalculate(p, p.x0)
	fmt.Println("Using Numerical Derivative Methods, solution  : ", solutionWithMethods)
	fmt.Println("Real solution : ", realSolution)
}
func GaussianJordanElimination() { //Inverse Matrix using GaussianJordanElimination
	p := Matrix{}
	p.takeMatrixSize()
	p.takeMatrixElements()
	p.displayMatrix()
	deger := p.determinat()
	fmt.Println(deger)
	if p.determinat() != 0 {
		unitMatrix := [10][10]float32{}
		//Function Closures
		unitMatrixConf := func() {
			for i := 0; i < 10; i++ {
				unitMatrix[i][i] = 1
			}
		}
		displayUnitMatrix := func() {
			fmt.Println("*** Display Matrix Inverse ***")
			for i := 0; i < p.matrixSize; i++ {
				for j := 0; j < p.matrixSize; j++ {
					fmt.Printf("%f  ", unitMatrix[i][j])
				}
				fmt.Println()
			}
		}
		unitMatrixConf()
		displayUnitMatrix()
		var tmp float32
		for i := 0; i < p.matrixSize; i++ {
			tmp = p.matrixElements[i][i]
			if tmp != 1.0 {
				for j := 0; j < p.matrixSize; j++ {
					p.matrixElements[i][j] /= tmp
					unitMatrix[i][j] /= tmp
				}
			}
			for j := 0; j < p.matrixSize; j++ {
				if i != j {
					tmp = p.matrixElements[j][i]
					for k := 0; k < p.matrixSize; k++ {
						p.matrixElements[j][k] -= tmp * p.matrixElements[i][k]
						unitMatrix[j][k] -= tmp * unitMatrix[i][k]
					}
				}
			}
		}
		displayUnitMatrix()
	} else {
		fmt.Println("None inverse of matrix,because determinat is zero(0) .")
	}
}
func GaussianElimination() {
	p := Matrix{}
	p.takeMatrixSize()
	p.takeMatrixElements()
	p.takeArrayElements()
	p.displayMatrix()
	p.displayArray()
	if p.determinat() != 0 { //Controlling Determinat
		var solution = [10]float32{1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
		//Function Closure
		displaySolution := func() {
			fmt.Println("*** Solution ***")
			fmt.Println("x1  x2  x3  x4  ....")
			for i := 0; i < p.matrixSize; i++ {
				fmt.Printf("%f  ", solution[i])
			}
			fmt.Println()
		}
		var tmp float32
		for i := 0; i < p.matrixSize; i++ {
			tmp = p.matrixElements[i][i]
			if tmp != 1.0 { //Controlling
				for j := 0; j < p.matrixSize; j++ {
					p.matrixElements[i][j] /= tmp
					p.arrayElements[i] /= tmp
				}
			}
			for j := i + 1; j < p.matrixSize; j++ {
				tmp = p.matrixElements[j][i]
				for k := i; k < p.matrixSize; k++ {
					p.matrixElements[j][k] /= tmp
					p.matrixElements[j][k] -= p.matrixElements[i][k]
					p.arrayElements[j] /= tmp
					p.arrayElements[j] -= p.arrayElements[i]
				}
			}
		}
		p.displayMatrix()
		p.displayArray()
		for i := p.matrixSize - 1; i >= 0; i-- {
			tmp = 0
			for j := 0; j < p.matrixSize; j++ {
				if i != j {
					tmp += solution[j] * p.matrixElements[i][j]
				}
			}
			solution[i] = p.arrayElements[i] - tmp
		}
		displaySolution()
	} else {
		fmt.Println("None inverse of matrix,because determinat is zero(0) .")
	}
}
func Trapez() {
	p := Necessary{}
	p.takeSize()
	p.takeElements()
	fmt.Println("First Range is top of integral, Second Range is bottom of integral")
	p.takeRanges()
	var numberSplit float32
	var solution float32
	fmt.Println("The number of split :")
	fmt.Scan(&numberSplit)
	p.deltaX = (p.firstRange - p.secondRange) / numberSplit
	for i := p.secondRange + p.deltaX; i < p.firstRange; i += p.deltaX {
		solution += fxCalculate(p, i)
	}
	solution = p.deltaX * ((fxCalculate(p, p.firstRange) + fxCalculate(p, p.secondRange)/2) + solution)
	fmt.Println("Solution : ", solution)

}
func Simpson() {
	p := Necessary{}
	p.takeSize()
	p.takeElements()
	fmt.Println("First Range is top of integral, Second Range is bottom of integral")
	p.takeRanges()
	var numberSplit float32
	var solution float32
	fmt.Println("The number of split :")
	fmt.Scan(&numberSplit)
	p.deltaX = (p.firstRange - p.secondRange) / numberSplit
	for i := p.secondRange + p.deltaX; i < p.firstRange; i += p.deltaX * 2 {
		solution += 4 * fxCalculate(p, i)
	}
	for i := p.secondRange + 2*p.deltaX; i < p.firstRange; i += p.deltaX * 2 {
		solution += 2 * fxCalculate(p, i)
	}
	solution = (p.deltaX / 3) * (fxCalculate(p, p.firstRange) + fxCalculate(p, p.secondRange) + solution)
	fmt.Println("Solution : ", solution)
}
func Secant() {
	p := Necessary{}
	p.takeSize()
	p.takeElements()
	p.takeX0()
	p.takeX1()
	p.takeEpsilon()
	p.x2 = p.x0 - (p.x1-p.x0)*fxCalculate(p, p.x0)/(fxCalculate(p, p.x1)-fxCalculate(p, p.x0))
	for Abs(p.x2-p.x1) > p.epsilon {
		p.x0 = p.x2
		p.x2 = p.x0 - (p.x1-p.x0)*fxCalculate(p, p.x0)/(fxCalculate(p, p.x1)-fxCalculate(p, p.x0))
	}
	fmt.Println("Root of equation : ", p.x2)
}
func NewtonRaphson() {
	p := Necessary{}
	p.takeSize()
	p.takeElements()
	p.takeX0()
	p.takeEpsilon()
	p.x0new = p.x0 - (fxCalculate(p, p.x0) / derivativeCalculate(p, p.x0))
	for Abs(p.x0new-p.x0) >= p.epsilon {
		p.x0 = p.x0new
		p.x0new = p.x0 - (fxCalculate(p, p.x0) / derivativeCalculate(p, p.x0))
	}
	fmt.Println("Root of equation : ", p.x0new)
}
func RegulaFalsi() {
	p := Necessary{}
	p.takeSize()
	p.takeElements()
	p.takeEpsilon()
	p.takeRanges()
	root1, root2 := fxCalculate(p, p.firstRange), fxCalculate(p, p.secondRange)
	for root1*root2 > 0 {
		p.takeRanges()
		root1, root2 = fxCalculate(p, p.firstRange), fxCalculate(p, p.secondRange)
	}
	fmt.Printf("\ta\t\tb\t   x=(b*f(a)-a*f(a))/(f(a)-f(b))\tf(a)\t\tf(b)\t\tf(x)\n\n")
	if root3 := p.epsilon + 1; root1*root2 != 0 {
		for Abs(root3) > p.epsilon {
			root1, root2 = fxCalculate(p, p.firstRange), fxCalculate(p, p.secondRange)
			p.thirdRange = (p.secondRange*root1 - p.firstRange*root2) / (root1 - root2)
			root3 = fxCalculate(p, p.thirdRange)
			fmt.Printf("%.8f  ||    %.8f   ||   %.25f    ||  %.8f   || %.8f ||  %.8f  ", p.firstRange, p.secondRange, p.thirdRange, root1, root2, root3)
			if root1*root3 < 0 {
				p.secondRange = p.thirdRange
			} else {
				p.firstRange = p.thirdRange
			}
			fmt.Println()
		}
		fmt.Println("Root of equation: ", p.thirdRange)
	} else {
		fmt.Println("Your Ranges are root or One of both is root.")
	}
}
func BiSection() {
	p := Necessary{}
	p.takeSize()
	p.takeElements()
	p.takeEpsilon()
	p.takeRanges()
	root1, root2 := fxCalculate(p, p.firstRange), fxCalculate(p, p.secondRange)
	for root1*root2 > 0 {
		p.takeRanges()
		root1, root2 = fxCalculate(p, p.firstRange), fxCalculate(p, p.secondRange)
	}
	fmt.Printf("\ta\t\tb\t\tx=(a+b)/2\tf(a)\t\tf(b)\t\tf(x)\n\n")
	if root3 := p.epsilon + 1; root1*root2 != 0 {
		for Abs(root3) > p.epsilon {
			root1, root2 = fxCalculate(p, p.firstRange), fxCalculate(p, p.secondRange)
			p.thirdRange = (p.firstRange + p.secondRange) / 2
			root3 = fxCalculate(p, p.thirdRange)
			fmt.Printf("%.8f  ||    %.8f   ||   %.8f    ||  %.8f   || %.8f ||  %.8f  ", p.firstRange, p.secondRange, p.thirdRange, root1, root2, root3)
			if root1*root3 < 0 {
				p.secondRange = p.thirdRange
			} else {
				p.firstRange = p.thirdRange
			}
			fmt.Println()
		}
		fmt.Println("Root of equation: ", p.thirdRange)
	} else {
		fmt.Println("Your Ranges are root or One of both is root.")
	}
}
func GraphicalAnalysis() { //Graphical Analysis Functions to finding root for equation
	r := Necessary{}
	r.takeSize()
	r.takeElements()
	r.takeX0()
	r.takeDeltaX()
	r.takeEpsilon()
	r.x0new = r.x0 + r.deltaX
	for (Abs(r.x0new-r.x0) >= r.epsilon) && fxCalculate(r, r.x0) != 0 {
		if fxCalculate(r, r.x0) < 0 {
			for fxCalculate(r, r.x0) < 0 {
				r.x0 += r.deltaX
			}
			r.x0new = r.x0
			r.x0 -= r.deltaX
			r.deltaX /= 2
		} else {
			for fxCalculate(r, r.x0) > 0 {
				r.x0 += r.deltaX
			}
			r.x0new = r.x0
			r.x0 -= r.deltaX
			r.deltaX /= 2
		}
	}
	fmt.Println("Root of equation : ", r.x0)
}

func Menu() int { //Display Menu
	fmt.Printf("\n\t\t<<<<MENU>>>>\n")
	fmt.Printf("\t[A] Root-Finding Methods\n")
	fmt.Println("[1]  Grafical")
	fmt.Println("[2]  Bi-Section ")
	fmt.Println("[3]  Regula-False ")
	fmt.Println("[4]  Newton-Raphson ")
	fmt.Println("[5]  Secant yontemi ")

	fmt.Printf("\t[B] Lineer Equation And Matrix\n")
	fmt.Println("[6]  Inverse Matrix with Gauss-Jordan elimination ")
	fmt.Println("[7]  Gauss elimination ")
	//	fmt.Println("[8]  Gauss Seidel")

	fmt.Printf("\t[C] Numerical Integral\n")
	fmt.Println("[9]  Trapez")
	fmt.Println("[10] Simpson")
	fmt.Printf("\t[D] Numerical Derivative Methods\n")
	fmt.Println("[11] Central Difference")
	fmt.Println("[0]  Exit")
	fmt.Print("Seciminizi giriniz:")
	var choice int
	fmt.Scanf("%d", &choice) //n, err := fmt.Scanln(&choice) //fmt.Scan(&choice)
	return choice
}
