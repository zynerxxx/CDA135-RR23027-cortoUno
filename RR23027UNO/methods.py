# 📚 Métodos para resolver sistemas de ecuaciones lineales y no lineales
# -------------------------------------------------------------
# Este archivo contiene la implementación de varios métodos matemáticos
# para resolver sistemas de ecuaciones. Cada método está documentado
# con su propósito, parámetros y valores de retorno.

import numpy as np

# ✨ Método: Eliminación de Gauss
# -------------------------------------------------------------
def eliminacion_gauss(A, b):
    """
    Resuelve el sistema lineal Ax = b utilizando el método de Eliminación de Gauss.
    
    Parámetros:
        A (lista de listas de floats): Matriz de coeficientes.
        b (lista de floats): Vector del lado derecho.

    Retorna:
        lista de floats: Vector solución x.
    """
    n = len(A)
    for i in range(n):
        # Pivot parcial
        max_row = max(range(i, n), key=lambda r: abs(A[r][i]))
        A[i], A[max_row] = A[max_row], A[i]
        b[i], b[max_row] = b[max_row], b[i]

        # Hacer el elemento diagonal 1
        diag = A[i][i]
        A[i] = [a / diag for a in A[i]]
        b[i] /= diag

        # Eliminar la columna debajo del pivote
        for j in range(i + 1, n):
            factor = A[j][i]
            A[j] = [a - factor * ai for a, ai in zip(A[j], A[i])]
            b[j] -= factor * b[i]

    # Sustitución hacia atrás
    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = b[i] - sum(A[i][j] * x[j] for j in range(i + 1, n))

    # Convertir los resultados a tipos estándar de Python
    return [float(xi) for xi in x]

# ✨ Método: Gauss-Jordan
# -------------------------------------------------------------
def gauss_jordan(A, b):
    """
    Resuelve el sistema lineal Ax = b utilizando el método de Gauss-Jordan.
    
    Parámetros:
        A (lista de listas de floats): Matriz de coeficientes.
        b (lista de floats): Vector del lado derecho.

    Retorna:
        lista de floats: Vector solución x.
    """
    n = len(A)
    for i in range(n):
        # Hacer el elemento diagonal 1
        diag = A[i][i]
        A[i] = [a / diag for a in A[i]]
        b[i] /= diag

        # Eliminar la columna para todas las demás filas
        for j in range(n):
            if i != j:
                factor = A[j][i]
                A[j] = [a - factor * ai for a, ai in zip(A[j], A[i])]
                b[j] -= factor * b[i]

    # Convertir los resultados a tipos estándar de Python
    return [float(xi) for xi in b]

# ✨ Método: Regla de Cramer
# -------------------------------------------------------------
def regla_cramer(A, b):
    """
    Resuelve el sistema lineal Ax = b utilizando la Regla de Cramer.
    
    Parámetros:
        A (lista de listas de floats): Matriz de coeficientes.
        b (lista de floats): Vector del lado derecho.

    Retorna:
        lista de floats: Vector solución x.
    """
    det_A = np.linalg.det(A)
    if det_A == 0:
        raise ValueError("El sistema no tiene una solución única.")

    n = len(A)
    x = []
    for i in range(n):
        A_copy = [row[:] for row in A]
        for row_index in range(n):
            A_copy[row_index][i] = b[row_index]
        x.append(np.linalg.det(A_copy) / det_A)

    # Convertir los resultados a tipos estándar de Python
    return [float(xi) for xi in x]

# ✨ Método: Descomposición LU
# -------------------------------------------------------------
def descomposicion_lu(A):
    """
    Realiza la Descomposición LU de una matriz A.

    Parámetros:
        A (lista de listas de floats): Matriz de coeficientes.

    Retorna:
        tupla: (L, U) donde L es la matriz triangular inferior y U es la matriz triangular superior.
    """
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]

    for i in range(n):
        # Triangular superior
        for k in range(i, n):
            U[i][k] = A[i][k] - sum(L[i][j] * U[j][k] for j in range(i))

        # Triangular inferior
        for k in range(i, n):
            if i == k:
                L[i][i] = 1  # Diagonal como 1
            else:
                L[k][i] = (A[k][i] - sum(L[k][j] * U[j][i] for j in range(i))) / U[i][i]

    return L, U

# ✨ Método: Jacobi
# -------------------------------------------------------------
def metodo_jacobi(A, b, x0, tol=1e-10, max_iterations=100):
    """
    Resuelve el sistema lineal Ax = b utilizando el método iterativo de Jacobi.

    Parámetros:
        A (lista de listas de floats): Matriz de coeficientes.
        b (lista de floats): Vector del lado derecho.
        x0 (lista de floats): Suposición inicial para la solución.
        tol (float): Tolerancia para la convergencia.
        max_iterations (int): Número máximo de iteraciones.

    Retorna:
        lista de floats: Vector solución x.
    """
    n = len(A)
    x = x0[:]
    for _ in range(max_iterations):
        x_new = [0] * n
        for i in range(n):
            s = sum(A[i][j] * x[j] for j in range(n) if j != i)
            x_new[i] = (b[i] - s) / A[i][i]
        if all(abs(x_new[i] - x[i]) < tol for i in range(n)):
            # Convertir los resultados a tipos estándar de Python
            return [float(xi) for xi in x_new]
        x = x_new
    raise ValueError("El método de Jacobi no convergió")

# ✨ Método: Gauss-Seidel
# -------------------------------------------------------------
def metodo_gauss_seidel(A, b, x0, tol=1e-10, max_iterations=100):
    """
    Resuelve el sistema lineal Ax = b utilizando el método iterativo de Gauss-Seidel.

    Parámetros:
        A (lista de listas de floats): Matriz de coeficientes.
        b (lista de floats): Vector del lado derecho.
        x0 (lista de floats): Suposición inicial para la solución.
        tol (float): Tolerancia para la convergencia.
        max_iterations (int): Número máximo de iteraciones.

    Retorna:
        lista de floats: Vector solución x.
    """
    n = len(A)
    x = x0[:]
    for _ in range(max_iterations):
        x_new = x[:]
        for i in range(n):
            s1 = sum(A[i][j] * x_new[j] for j in range(i))
            s2 = sum(A[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (b[i] - s1 - s2) / A[i][i]
        if all(abs(x_new[i] - x[i]) < tol for i in range(n)):
            # Convertir los resultados a tipos estándar de Python
            return [float(xi) for xi in x_new]
        x = x_new
    raise ValueError("El método de Gauss-Seidel no convergió")

# ✨ Método: Bisección
# -------------------------------------------------------------
def metodo_biseccion(f, a, b, tol=1e-10, max_iterations=100):
    """
    Encuentra una raíz de la ecuación f(x) = 0 utilizando el método de Bisección.

    Parámetros:
        f (función): La función para la cual se busca la raíz.
        a (float): El límite inferior del intervalo.
        b (float): El límite superior del intervalo.
        tol (float): Tolerancia para la convergencia.
        max_iterations (int): Número máximo de iteraciones.

    Retorna:
        float: La raíz de la función.
    """
    if f(a) * f(b) >= 0:
        raise ValueError("La función debe tener signos opuestos en los extremos a y b.")

    for _ in range(max_iterations):
        c = (a + b) / 2
        if abs(f(c)) < tol or (b - a) / 2 < tol:
            return float(c)
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    raise ValueError("El método de Bisección no convergió")