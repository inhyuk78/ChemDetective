# Handling file export 
def export_to_csv(df):
    df.to_csv('output_file.csv', index=True)
    print('Results exported as output_file.csv file in the directory.')

def export_to_excel(df):
    df.to_excel('output_file.xlsx')
    print('Results exported as output_file.xlsx file in the directory.')
