import math

a = 6378137.0  # Earth's equatorial radius in meters
f = 1 / 298.257223563  # Earth's flattening factor


def vincenty_distance(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """
    Calculate the geodesic distance between two points on the Earth's surface using the Vincenty formula.
    
    The Vincenty formula is an accurate method to calculate the shortest distance between two points
    on the surface of an ellipsoid, such as the Earth. It accounts for the ellipsoidal nature of the Earth,
    providing higher accuracy compared to simpler formulas like the Haversine formula. The formula involves
    an iterative process to refine the distance calculation.
    
    Args:
        lat1 (float): Latitude of the first point in degrees.
        lon1 (float): Longitude of the first point in degrees.
        lat2 (float): Latitude of the second point in degrees.
        lon2 (float): Longitude of the second point in degrees.
    
    Returns:
        float: Geodesic distance between the two points in meters.
    
    Example:
        >>> lat1 = 34.0522  # Los Angeles latitude in degrees
        >>> lon1 = -118.2437  # Los Angeles longitude in degrees
        >>> lat2 = 40.7128  # New York City latitude in degrees
        >>> lon2 = -74.0060  # New York City longitude in degrees
        >>> distance = vincenty_distance(lat1, lon1, lat2, lon2)
        >>> print(f"Distance between Los Angeles and New York City: {distance:.2f} meters")
    """
  
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    lambda1 = math.radians(lon1)
    lambda2 = math.radians(lon2)
    
    U1 = math.atan((1 - f) * math.tan(phi1))
    U2 = math.atan((1 - f) * math.tan(phi2))
    L = lambda2 - lambda1
    lambda_ = L
    
    iter_limit = 1000
    for _ in range(iter_limit):
        sin_lambda = math.sin(lambda_)
        cos_lambda = math.cos(lambda_)
        sin_sigma = math.sqrt((math.cos(U2) * sin_lambda) ** 2 +
                              (math.cos(U1) * math.sin(U2) -
                               math.sin(U1) * math.cos(U2) * cos_lambda) ** 2)
        cos_sigma = math.sin(U1) * math.sin(U2) + math.cos(U1) * math.cos(U2) * cos_lambda
        sigma = math.atan2(sin_sigma, cos_sigma)
        sin_alpha = math.cos(U1) * math.cos(U2) * sin_lambda / sin_sigma
        cos2_alpha = 1 - sin_alpha ** 2
        cos2_sigma_m = cos_sigma - 2 * math.sin(U1) * math.sin(U2) / cos2_alpha
        C = f / 16 * cos2_alpha * (4 + f * (4 - 3 * cos2_alpha))
        lambda_prev = lambda_
        lambda_ = L + (1 - C) * f * sin_alpha * (
                sigma + C * sin_sigma * (
                cos2_sigma_m + C * cos_sigma * (-1 + 2 * cos2_sigma_m ** 2)))
        if abs(lambda_ - lambda_prev) < 1e-12:
            break
    
    u2 = cos2_alpha * (a ** 2 - b ** 2) / (b ** 2)
    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    delta_sigma = B * sin_sigma * (
            cos2_sigma_m + B / 4 * (
            cos_sigma * (-1 + 2 * cos2_sigma_m ** 2) - B / 6 * cos2_sigma_m * (
            -3 + 4 * sin_sigma ** 2) * (-3 + 4 * cos2_sigma_m ** 2)))
    distance = b * A * (sigma - delta_sigma)
    
    return distance

if __name__=='__main__':
  lat1 = 34.0522  # Los Angeles latitude in degrees
  lon1 = -118.2437  # Los Angeles longitude in degrees
  lat2 = 40.7128  # New York City latitude in degrees
  lon2 = -74.0060  # New York City longitude in degrees

  distance = vincenty_distance(lat1, lon1, lat2, lon2)
  print(f"Distance between Los Angeles and New York City: {distance:.2f} meters")
