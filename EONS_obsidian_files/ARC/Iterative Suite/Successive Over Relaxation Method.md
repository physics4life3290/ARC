
The **Successive Over-Relaxation (SOR)** method is an iterative technique used for solving **linear systems of equations**. It is a refinement of the **Gauss-Seidel method**, aimed at improving convergence speed by introducing a relaxation factor.

The basic idea of the **SOR method** is to adjust the Gauss-Seidel iteration by a factor \( \omega \) (known as the **relaxation parameter**). This parameter can accelerate the convergence of the method by allowing more or less "weight" to be placed on the most recent updates to the solution vector.

## **1. Overview of the SOR Method**

The SOR method solves a system of linear equations \( Ax = b \), where \( A \) is a matrix, \( x \) is the vector of unknowns, and \( b \) is the right-hand side vector.

Given a linear system:

$$
A x = b
$$

The SOR method refines the **Gauss-Seidel method** by introducing the **relaxation parameter** \( \omega \), where \( 0 < \omega < 2 \). The update rule for the solution vector \( x \) is:

$$
x_i^{(k+1)} = (1 - \omega) x_i^{(k)} + \frac{\omega}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij} x_j^{(k)} \right)
$$

Where:
- \( x_i^{(k+1)} \) is the updated value of the \( i \)-th unknown at the \( (k+1) \)-th iteration,
- \( x_i^{(k)} \) is the current value of the \( i \)-th unknown at the \( k \)-th iteration,
- \( a_{ii} \) is the diagonal entry of matrix \( A \) for the \( i \)-th equation,
- \( b_i \) is the right-hand side of the equation.

The key feature of the SOR method is the **relaxation factor** \( \omega \), which modifies the Gauss-Seidel updates to accelerate convergence.

## **2. How the SOR Method Works**

1. **Initialization**: Start with an initial guess for the solution vector \( x^{(0)} \).
2. **Iteration**: For each iteration \( k \), update each component of \( x_i \) using the formula:
   $$
   x_i^{(k+1)} = (1 - \omega) x_i^{(k)} + \omega \cdot \left( \text{Gauss-Seidel update for } x_i \right)
   $$ 
   The update for \( x_i \) incorporates the previously computed values in the current iteration, but the relaxation factor \( \omega \) controls how much weight is placed on the new value versus the old one.
3. **Convergence Check**: After each iteration, check for convergence by comparing the difference between the new and old solutions. If the solution has converged within a specified tolerance, stop the iteration.

## **3. Choosing the Relaxation Parameter \( \omega \)**

The convergence of the **SOR method** depends heavily on the value of \( \omega \). If \( \omega = 1 \), the SOR method reduces to the **Gauss-Seidel method**. By choosing values of \( \omega \) between 0 and 2, the convergence can be accelerated, but it is crucial to find an optimal \( \omega \) for the specific problem.

- If \( \omega \) is too small, the method converges slowly (similar to Gauss-Seidel).
- If \( \omega \) is too large, the method may **diverge** or oscillate around the solution.

Optimal values of \( \omega \) are often found experimentally or by using **spectral analysis**. Typical values of \( \omega \) for practical problems range from 1.0 to 1.5.

## **4. Convergence Criteria**

For the SOR method to converge, the matrix \( A \) needs to satisfy certain properties:

1. **Symmetric Positive-Definite Matrix**: If the matrix \( A \) is symmetric and positive-definite, the SOR method will converge for any initial guess, given an appropriate value of \( \omega \).
2. **Diagonally Dominant Matrix**: For matrices that are not symmetric, the SOR method may still converge if the matrix is **diagonally dominant** (i.e., the magnitude of each diagonal entry is greater than the sum of the magnitudes of the other entries in the corresponding row).

The convergence condition can be written as:

$$
|a_{ii}| > \sum_{j \neq i} |a_{ij}|
$$

where the matrix \( A \) is diagonally dominant.

## **5. Advantages of the SOR Method**

### **5.1. Faster Convergence**

The **SOR method** can converge faster than the **Gauss-Seidel method**, especially when a well-chosen relaxation parameter \( \omega \) is used. This makes it particularly useful for large systems of equations, such as those arising in **finite element analysis** or **computational fluid dynamics**.

### **5.2. Flexible and Parallelizable**

Since the SOR method involves updating each unknown sequentially, it is highly parallelizable, making it well-suited for modern **parallel computing** environments. Additionally, the method is flexible and can be adapted to a variety of problem types.

### **5.3. Memory Efficient**

The **SOR method** requires only the storage of the solution vector and matrix coefficients, making it memory efficient compared to direct methods such as **Gaussian elimination**, which require storing the entire augmented matrix.

## **6. Disadvantages of the SOR Method**

### **6.1. Sensitive to Initial Guess**

The convergence of the **SOR method** depends on the initial guess \( x^{(0)} \). If the initial guess is far from the true solution, the method may converge slowly or fail to converge. Proper preconditioning or a good initial guess is crucial for efficient convergence.

### **6.2. Selection of Relaxation Parameter**

Choosing the **optimal relaxation parameter \( \omega \)** can be difficult. A poor choice can lead to slow convergence or even divergence. In practice, a range of \( \omega \) values is often tested to find the one that minimizes the number of iterations.

### **6.3. Not Suitable for All Problems**

The **SOR method** is not suitable for all types of systems. For example, it may not be effective for highly ill-conditioned matrices or systems with complex boundary conditions.

## **7. Applications of the SOR Method**

The **SOR method** is widely used in fields such as:

- **Computational Fluid Dynamics (CFD)**: For solving large systems of linear equations that arise in fluid flow simulations.
- **Finite Element Analysis (FEA)**: Used to solve large sparse systems of linear equations in structural engineering problems.
- **Image Reconstruction**: Applied in image processing for tasks such as medical imaging or signal processing.
- **Astrophysical Simulations**: Used in solving linear systems that arise in large-scale simulations, such as **core-collapse supernova models** or **stellar evolution models**.

## **8. Summary**

- The **Successive Over-Relaxation (SOR) method** accelerates the **Gauss-Seidel method** by introducing a relaxation parameter \( \omega \).
- The method updates each variable using the weighted average of the old and new values, where the relaxation parameter controls the update.
- The method requires careful selection of \( \omega \) and is typically most effective for **symmetric positive-definite** or **diagonally dominant** systems.
- It is widely used in **numerical simulations**, particularly for solving large systems in **fluid dynamics**, **finite element analysis**, and **astrophysical models**.

---

Tags:
[[Iterative Suite]]