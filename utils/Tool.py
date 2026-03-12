import ROOT
import sys

def get_sparse(file_list, name):
        
        sparse = []
        for file_dir in file_list:
            file = ROOT.TFile(file_dir, "READ")
            sparse.append(file.Get(name))
            file.Close()
      
        for i in range(1, len(sparse)):
            try:
                sparse[0].Add(sparse[i])
                sparse[i].Delete()
            except:
                print(f"Error adding sparse from file: {file_list[i]}")
                sys.exit(1)
            
        return sparse[0]