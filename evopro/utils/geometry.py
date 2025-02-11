import math

def distance(p1, p2):
    """returns the distance between two 3D points represented as tuples"""

    dist = math.sqrt((float(p2[0])-float(p1[0]))**2+(float(p2[1])-(float(p1[1])))**2+(float(p2[2])-float(p1[2]))**2)
    return dist