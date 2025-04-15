# 📘 Biblioteca RR23027UNO

## ✨ Descripción
La biblioteca `RR23027UNO` proporciona métodos para resolver sistemas de ecuaciones lineales y no lineales. Incluye los siguientes métodos:

- 🔹 **Eliminación de Gauss**
- 🔹 **Gauss-Jordan**
- 🔹 **Método de Cramer**
- 🔹 **Descomposición LU**
- 🔹 **Método de Jacobi**
- 🔹 **Método de Gauss-Seidel**
- 🔹 **Método de Bisección**

Esta biblioteca está diseñada para ser fácil de usar y comprender, se trabaja con métodos numéricos.

---

## 🛠️ Instalación
Para instalar la biblioteca, utilizaremos el siguiente comando:

```bash
pip install RR23027UNO
```

---

## 📋 Uso
Aquí comprendemos un ejemplo de cómo usar la biblioteca:

```python
from RR23027UNO.methods import eliminacion_gauss

A = [[2, -1, 1], [3, 3, 9], [3, 3, 5]]
b = [2, -1, 4]
solucion = eliminacion_gauss(A, b)
print("Solución:", solucion)
```

Para más ejemplos, podemos consultar el archivo `ejemploUso.py` (con py).

---

## 📂 Estructura del Proyecto
```
RR23027UNO/
├── __init__.py
├── methods.py
├── __pycache__/
└── ejemploUso.py
```

---

## 💡 Notas
- Asegurarnos de tener instalada la biblioteca `numpy` antes de usar los métodos.
