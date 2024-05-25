CREATE TABLE Alumno (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    apellido VARCHAR(255),
    -- otros atributos
);

CREATE TABLE Libro (
    id INT PRIMARY KEY,
    titulo VARCHAR(255),
    autor VARCHAR(255),
    -- otros atributos
);

CREATE TABLE Prestamo (
    id INT PRIMARY KEY,
    libro_id INT,
    alumno_id INT,
    fecha_prestamo DATE,
    fecha_devolucion DATE,
    FOREIGN KEY (libro_id) REFERENCES Libro(id),
    FOREIGN KEY (alumno_id) REFERENCES Alumno(id)
);

CREATE TABLE Revista (
    id INT PRIMARY KEY,
    titulo VARCHAR(255),
    -- otros atributos
);

CREATE TABLE Articulo (
    id INT PRIMARY KEY,
    titulo VARCHAR(255),
    revista_id INT,
    FOREIGN KEY (revista_id) REFERENCES Revista(id)
);

CREATE TABLE Sucursal (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    -- otros atributos
);

CREATE TABLE Vino (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    -- otros atributos
);

CREATE TABLE Componente (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    porcentaje FLOAT,
    vino_id INT,
    FOREIGN KEY (vino_id) REFERENCES Vino(id)
);

CREATE TABLE Cliente (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    -- otros atributos
);

CREATE TABLE Accidente (
    id INT PRIMARY KEY,
    cliente_id INT,
    fecha DATE,
    FOREIGN KEY (cliente_id) REFERENCES Cliente(id)
);

CREATE TABLE Practica (
    id INT PRIMARY KEY,
    descripcion VARCHAR(255),
    -- otros atributos
);

CREATE TABLE Doctor (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    apellido VARCHAR(255),
    -- otros atributos
);

CREATE TABLE PracticaRealizada (
    id INT PRIMARY KEY,
    practica_id INT,
    doctor_id INT,
    paciente_id INT,
    fecha DATE,
    FOREIGN KEY (practica_id) REFERENCES Practica(id),
    FOREIGN KEY (doctor_id) REFERENCES Doctor(id),
    FOREIGN KEY (paciente_id) REFERENCES Cliente(id) -- Asumiendo que los pacientes son clientes
);

CREATE TABLE Proveedor (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    -- otros atributos
);

CREATE TABLE Ingrediente (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    -- otros atributos
);

CREATE TABLE ProveedorIngrediente (
    id INT PRIMARY KEY,
    proveedor_id INT,
    ingrediente_id INT,
    cantidad INT,
    FOREIGN KEY (proveedor_id) REFERENCES Proveedor(id),
    FOREIGN KEY (ingrediente_id) REFERENCES Ingrediente(id)
);

CREATE TABLE Producto (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    -- otros atributos
);

CREATE TABLE IngredienteProducto (
    id INT PRIMARY KEY,
    ingrediente_id INT,
    producto_id INT,
    cantidad INT,
    FOREIGN KEY (ingrediente_id) REFERENCES Ingrediente(id),
    FOREIGN KEY (producto_id) REFERENCES Producto(id)
);

CREATE TABLE Proyecto (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    presupuesto DECIMAL(10, 2),
    gerencia_id INT,
    FOREIGN KEY (gerencia_id) REFERENCES Gerencia(id)
);

CREATE TABLE Gerencia (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    -- otros atributos
);

CREATE TABLE Guia (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    apellido VARCHAR(255),
    -- otros atributos
);

CREATE TABLE Tour (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    guia_id INT,
    FOREIGN KEY (guia_id) REFERENCES Guia(id)
);

CREATE TABLE Pasajero (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    apellido VARCHAR(255),
    -- otros atributos
);

CREATE TABLE TourPasajero (
    id INT PRIMARY KEY,
    tour_id INT,
    pasajero_id INT,
    precio DECIMAL(10, 2),
    FOREIGN KEY (tour_id) REFERENCES Tour(id),
    FOREIGN KEY (pasajero_id) REFERENCES Pasajero(id)
);

CREATE TABLE Equipo (
    id INT PRIMARY KEY,
    nombre VARCHAR(255),
    -- otros atributos
);

CREATE TABLE Partido (
    id INT PRIMARY KEY,
    equipo_local_id INT,
    equipo_visitante_id INT,
    goles_local INT,
    goles_visitante INT,
    fecha DATE,
    FOREIGN KEY (equipo_local_id) REFERENCES Equipo(id),
    FOREIGN KEY (equipo_visitante_id) REFERENCES Equipo(id)
);


-- Consulta de préstamos que ya deberían estar devueltos
SELECT l.titulo AS libro, a.nombre AS alumno, a.apellido AS alumno_apellido, p.fecha_prestamo, p.fecha_devolucion
FROM Prestamo p
JOIN Libro l ON p.libro_id = l.id
JOIN Alumno a ON p.alumno_id = a.id
WHERE p.fecha_devolucion < CURRENT_DATE AND p.fecha_devolucion IS NOT NULL;

-- Consulta de la cantidad de artículos por ejemplar de la revista "National Geographics"
SELECT r.titulo AS revista, COUNT(a.id) AS cantidad_articulos
FROM Revista r
JOIN Articulo a ON r.id = a.revista_id
WHERE r.titulo = 'National Geographics'
GROUP BY r.titulo;

-- Consulta de artículos vendidos por la sucursal 1 y no en la sucursal 2
SELECT a1.id, a1.titulo
FROM Articulo a1
JOIN SucursalArticulo sa1 ON a1.id = sa1.articulo_id
LEFT JOIN SucursalArticulo sa2 ON a1.id = sa2.articulo_id AND sa2.sucursal_id = 2
WHERE sa1.sucursal_id = 1 AND sa2.articulo_id IS NULL;

-- Consulta de vinos donde la sumatoria de % de componentes no es 100%
SELECT v.nombre, SUM(c.porcentaje) AS suma_porcentajes
FROM Vino v
JOIN Componente c ON v.id = c.vino_id
GROUP BY v.nombre
HAVING SUM(c.porcentaje) <> 100;

-- Consulta de clientes que nunca han tenido un accidente
SELECT c.id, c.nombre, c.apellido
FROM Cliente c
LEFT JOIN Accidente a ON c.id = a.cliente_id
WHERE a.id IS NULL;

-- Consulta de prácticas que el Dr. Cureta no ha realizado a ningún paciente
SELECT p.id, p.descripcion
FROM Practica p
LEFT JOIN PracticaRealizada pr ON p.id = pr.practica_id AND pr.doctor_id = (SELECT id FROM Doctor WHERE nombre = 'Cureta')
WHERE pr.id IS NULL;

-- Consulta del proveedor más importante (el que provee la mayor cantidad de ingredientes)
SELECT p.id, p.nombre, SUM(pi.cantidad) AS total_ingredientes
FROM Proveedor p
JOIN ProveedorIngrediente pi ON p.id = pi.proveedor_id
GROUP BY p.id, p.nombre
ORDER BY total_ingredientes DESC
LIMIT 1;

-- Consulta del total de $ involucrados en los proyectos de la gerencia de Marketing
SELECT SUM(pr.presupuesto) AS total_presupuesto
FROM Proyecto pr
JOIN Gerencia g ON pr.gerencia_id = g.id
WHERE g.nombre = 'Marketing';

-- Consulta del guía más significativo en cuanto a $ de ingresos
SELECT g.id, g.nombre, g.apellido, SUM(tp.precio) AS total_ingresos
FROM Guia g
JOIN Tour t ON g.id = t.guia_id
JOIN TourPasajero tp ON t.id = tp.tour_id
GROUP BY g.id, g.nombre, g.apellido
ORDER BY total_ingresos DESC
LIMIT 1;

-- Consulta del listado de partidos jugados con resultados
SELECT p.id, e_local.nombre AS equipo_local, e_visitante.nombre AS equipo_visitante,
       p.goles_local, p.goles_visitante,
       CASE
           WHEN p.goles_local > p.goles_visitante THEN 'Local'
           WHEN p.goles_local < p.goles_visitante THEN 'Visitante'
           ELSE 'Empate'
       END AS resultado
FROM Partido p
JOIN Equipo e_local ON p.equipo_local_id = e_local.id
JOIN Equipo e_visitante ON p.equipo_visitante_id = e_visitante.id;
