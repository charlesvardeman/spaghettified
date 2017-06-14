"""
.. module:: outputLogger
    :synopsis: Allows logging to the specified location and can overwrite or append as well as print to stdout.
    

.. moduleauthor:: Nathan Smith <nsmith13@nd.edu>

"""

class outputLogger ():
    def __init__(self,filePath=None,mode="overwrite",store=True,printout=True):
        """Initializes the logger.
        
        Kwargs:
            | filePath (str): 'dir/file.txt'
            | mode (str): 'overwrite' or 'append'
            | store (bool): Whether to store the messages in a file.
            | printout (bool): Whether to print the message to stdout.
        """
        
        self.store = store
        self.printout = printout
        self.filePath = filePath
        if self.filePath == None:
            self.filePath = 'outputs/logger_output.txt'
        if self.store == True and mode == "overwrite":
            f = file(self.filePath,'w')
            f.write('')
            f.close()
    
    def log(self, message):
        """Store and/or prints the message depending on options
        
        Args:
            message (object): Any object that has a __str__ method
        """
        
        if self.store:
            f = file(self.filePath,'a')
            f.write(str(message)+'\n')
            f.close()
        
        if self.printout:
            print message

