from RR23027UNO.methods import *

def main():
    print("\n==============================")
    print("✨ Ejemplo de Métodos Numéricos ✨")
    print("==============================\n")

    # Ejemplo para Eliminación de Gauss
    print("🔹 Eliminación de Gauss")
    A = [[2, -1, 1], [3, 3, 9], [3, 3, 5]]
    b = [2, -1, 4]
    solucion = eliminacion_gauss(A, b)
    print("Resultado:", solucion)
    print("------------------------------\n")

    # Ejemplo para Gauss-Jordan
    print("🔹 Gauss-Jordan")
    A = [[2, -1, 1], [3, 3, 9], [3, 3, 5]]
    b = [2, -1, 4]
    solucion = gauss_jordan(A, b)
    print("Resultado:", solucion)
    print("------------------------------\n")

    # Ejemplo para la Regla de Cramer
    print("🔹 Regla de Cramer")
    A = [[2, -1, 1], [3, 3, 9], [3, 3, 5]]
    b = [2, -1, 4]
    solucion = regla_cramer(A, b)
    print("Resultado:", solucion)
    print("------------------------------\n")

    # Ejemplo para Descomposición LU
    print("🔹 Descomposición LU")
    A = [[2, -1, 1], [3, 3, 9], [3, 3, 5]]
    L, U = descomposicion_lu(A)
    print("Matriz L:", L)
    print("Matriz U:", U)
    print("------------------------------\n")

    # Ejemplo para el Método de Jacobi
    print("🔹 Método de Jacobi")
    A = [[4, -1, 0, 0], [-1, 4, -1, 0], [0, -1, 4, -1], [0, 0, -1, 3]]
    b = [15, 10, 10, 10]
    x0 = [0, 0, 0, 0]
    solucion = metodo_jacobi(A, b, x0)
    print("Resultado:", solucion)
    print("------------------------------\n")

    # Ejemplo para el Método de Gauss-Seidel
    print("🔹 Método de Gauss-Seidel")
    solucion = metodo_gauss_seidel(A, b, x0)
    print("Resultado:", solucion)
    print("------------------------------\n")

    # Ejemplo para el Método de Bisección
    print("🔹 Método de Bisección")
    f = lambda x: x**3 - x - 2
    solucion = metodo_biseccion(f, 1, 2)
    print("Resultado:", solucion)
    print("==============================")

if __name__ == "__main__":
    main()