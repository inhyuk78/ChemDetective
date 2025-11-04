# Handling file export 
def export_to_csv(df):
    return df.to_csv('output_file.csv', index=True)

def export_to_excel(df):
    return df.to_excel('output_file.xlsx')
