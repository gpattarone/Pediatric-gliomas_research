#!/usr/bin/env nextflow

process LoadMetaAnalysis {
    input:
    path csv_file
    output:
    stdout result

    script:
    """
    python -c "
import pandas as pd
import psycopg2

def load_to_postgresql(csv_file):
    # Database connection details
    conn = psycopg2.connect(
        dbname='therapies_db',
        user='your_username',
        password='your_password',
        host='localhost',
        port='5432'
    )
    cursor = conn.cursor()

    # Define the schema
    create_table_query = '''
    CREATE TABLE IF NOT EXISTS meta_analysis (
        Therapy TEXT,
        Target_Group TEXT,
        Evidence TEXT,
        Comments TEXT
    )
    '''
    cursor.execute(create_table_query)

    # Load CSV into DataFrame
    df = pd.read_csv(csv_file)

    # Insert data into PostgreSQL
    for _, row in df.iterrows():
        insert_query = '''
        INSERT INTO meta_analysis (Therapy, Target_Group, Evidence, Comments)
        VALUES (%s, %s, %s, %s)
        '''
        cursor.execute(insert_query, (row['Therapy'], row['Target_Group'], row['Evidence'], row['Comments']))

    conn.commit()
    cursor.close()
    conn.close()

load_to_postgresql('${csv_file}')
"
    """
}

workflow {
    csv_input = file(params.meta_analysis_csv)
    LoadMetaAnalysis(csv_input)
}
