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

Esta biblioteca está diseñada para ser fácil de usar y comprender, ideal para estudiantes y profesionales que trabajan con métodos numéricos.

---

## 🛠️ Instalación
Para instalar la biblioteca, utiliza el siguiente comando después de publicarla en PyPI:

```bash
pip install RR23027UNO
```

---

## 📋 Uso
Aquí tienes un ejemplo de cómo usar la biblioteca:

```python
from RR23027UNO.methods import eliminacion_gauss

A = [[2, -1, 1], [3, 3, 9], [3, 3, 5]]
b = [2, -1, 4]
solucion = eliminacion_gauss(A, b)
print("Solución:", solucion)
```

Para más ejemplos, consulta el archivo `ejemploUso.py`.

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

## 📜 Licencia
Este proyecto está licenciado bajo la Licencia MIT.

---

## 💡 Notas
- Asegúrate de tener instalada la biblioteca `numpy` antes de usar los métodos.
- Si encuentras algún problema, no dudes en reportarlo.