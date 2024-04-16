import subprocess
import os

def descargar_archivos_sra(numero_acceso):
    print("Descargando archivos SRA...")
    subprocess.run(["prefetch", "-O", ".", numero_acceso])

def convertir_a_fastq(numero_acceso):
    print("Convirtiendo archivos SRA a formato FastQ...")
    subprocess.run(["fastq-dump", "--split-files", numero_acceso])

def analisis_con_rsem():
    print("Realizando análisis con RSEM...")
    subprocess.run(["rsem-calculate-expression", "-p", "8", "--paired-end", "sample_1.fastq", "sample_2.fastq", "reference_transcriptome", "output_prefix"])

def analisis_informacion_secuencias(archivo_fasta):
    print("Realizando análisis de información de secuencias...")
    # Aquí deberías escribir el código para analizar el archivo fasta y guardar la información relevante en transcript_information.txt

def analisis_expresion_transcriptos():
    print("Realizando análisis de expresión de transcriptos...")
    # Aquí deberías escribir el código para calcular los valores de RPKM y guardarlos en transcript_expression.txt

def main():
    # Funcionalidad 1: Descarga de Archivos SRA
    numero_acceso_sra = input("Ingrese el número de acceso SRA: ")
    descargar_archivos_sra(numero_acceso_sra)
    print("Archivos SRA descargados correctamente.")

    # Funcionalidad 2: Conversión a FastQ
    convertir_a_fastq(numero_acceso_sra)
    print("Archivos SRA convertidos a formato FastQ.")

    # Funcionalidad 3: Análisis con RSEM
    analisis_con_rsem()
    print("Análisis con RSEM completado.")

    # Funcionalidad 4: Análisis de Información de Secuencias
    archivo_fasta = input("Ingrese la ruta del archivo fasta: ")
    analisis_informacion_secuencias(archivo_fasta)
    print("Análisis de información de secuencias completado.")

    # Funcionalidad 5: Análisis de Expresión de Transcriptos
    analisis_expresion_transcriptos()
    print("Análisis de expresión de transcriptos completado.")

if __name__ == "__main__":
    main()
