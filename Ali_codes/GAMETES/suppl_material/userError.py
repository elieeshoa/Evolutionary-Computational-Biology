class userError(Exception):
    """ 
    A calss for returning user defined errors 

    Ali R. Zomorrodi, Segre lab @ Bosotn University
    """
    def __init__(self, error_msgs = ''):
        self.error = '\n**ERROR! ' + str(error_msgs) + '\n'

    def __str__(self):
        return self.error


