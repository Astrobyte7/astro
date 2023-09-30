from flask import Flask, render_template
from flask import Flask, render_template, request
from flask import Flask, render_template, request, jsonify
from astropy.coordinates import get_body_barycentric_posvel

import ephem
import datetime
from datetime import datetime
import pytz
from skyfield import api, almanac
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
from skyfield.api import Topos, load
import io
import base64
from skyfield import eclipselib
from skyfield.api import load
from astropy.time import Time
import geocoder
from geocoder import osm
from pytz import timezone
from skyfield.api import N, W, S, E, wgs84, load
import os
import csv
import folium
import math
from io import BytesIO

import requests


app = Flask(__name__)

@app.route('/')
def home():
    return render_template('home.html')

@app.route('/solar_system')
def solar_system():
    return render_template('solar_system.html')

@app.route('/constellations')
def constellations():
    return render_template('constellations.html')

@app.route('/blog')
def blog():
    return render_template('blog.html')

@app.route('/contact')
def contact():
    return render_template('contact.html')

@app.route('/About_me')
def About_me():
    return render_template('About_me.html')

@app.route('/astro_calender')
def astro_calender():
    return render_template('astro_calender.html')

@app.route('/learn')
def learn():
    return render_template('learn.html')



@app.route('/learn_solar')
def learn_solar():
    return render_template('learn_solar.html')

@app.route('/learn_astroid')
def learn_astroid():
    return render_template('learn_astroid.html')


@app.route('/learn_planets')
def learn_planets():
    return render_template('learn_planets.html')



@app.route('/learn_galaxy')
def learn_galaxy():
    return render_template('learn_galaxy.html')


@app.route('/blog_henlay')
def blog_henlay():
    return render_template('blog_henlay.html')


@app.route('/blog_carl')
def blog_carl():
    return render_template('blog_carl.html')




@app.route('/learn_astrobilogy')
def learn_astrobilogy():
    return render_template('learn_astrobilogy.html')


@app.route('/learn_supernova')
def learn_supernova():
    return render_template('learn_supernova.html')
















@app.route('/planets', methods=['GET', 'POST'])
def planets():
    results = None
    celestial_bodies = [
        'Sun',
        'Moon',
        'Mercury',
        'Venus',
        'Mars',
        'Jupiter',
        'Saturn',
        'Uranus',
        'Neptune',
    ]
    timezones = [
        ('Indian Standard Time', 'Asia/Kolkata'),
        ('Eastern Standard Time', 'America/New_York'),
        ('Coordinated Universal Time', 'UTC'),
        ('London Time', 'Europe/London'),
        # You can add more timezone options here
    ]

    date_str = dt.datetime.now().strftime('%Y-%m-%d')  # Use the current date by default
    timezone_str = 'Indian Standard Time'  # Set IST as default timezone
    city_str = 'Ahmedabad'  # Set Mumbai as default city

    if request.method == 'POST':
        date_str = request.form['date']
        timezone_str = request.form['timezone']
        city_str = request.form['city']

    date = datetime.strptime(date_str, '%Y-%m-%d')
    observer = ephem.Observer()

    location = geocoder.osm(city_str)
    if location.ok:
        observer.lat, observer.lon = str(location.latlng[0]), str(location.latlng[1])
    else:
        return render_template('error.html', message='Geocoding failed for the provided city.')

    for tz_name, tz in timezones:
        if tz_name == timezone_str:
            local_tz = pytz.timezone(tz)
            observer.date = local_tz.localize(date)
            break

    results = {}

    for celestial_body_name in celestial_bodies:
        celestial_body = getattr(ephem, celestial_body_name)()
        celestial_body.compute(observer)

        distance_sun = format(celestial_body.sun_distance, '.2f')
        distance_earth = format(celestial_body.earth_distance, '.2f')

        transit_time = observer.previous_transit(celestial_body).datetime().replace(tzinfo=pytz.UTC).astimezone(local_tz)

        results[celestial_body_name.lower()] = {
            'rise_time': observer.next_rising(celestial_body).datetime().replace(tzinfo=pytz.UTC).astimezone(local_tz),
            'set_time': observer.next_setting(celestial_body).datetime().replace(tzinfo=pytz.UTC).astimezone(local_tz),
            'transit_time': transit_time.strftime('%H:%M'),
            'phase': format(celestial_body.phase, '.2f'),
            'constellation': ephem.constellation(celestial_body),
            'magnitude': format(celestial_body.mag, '.2f'),
            'distance_sun': distance_sun,
            'distance_earth': distance_earth,
        }

    return render_template('planets.html', results=results, timezones=timezones, date_str=date_str, timezone_str=timezone_str, city_str=city_str)






def calculate_astronomical_data(date, city, selected_timezone):
    location = osm(city)
    if not location.ok:
        return None  # Handle the case when geocoding fails

    lat, lon = location.latlng

    home = ephem.Observer()
    home.date = date
    home.lat, home.lon = str(lat), str(lon)

    sun, moon = ephem.Sun(), ephem.Moon()

    sun.compute(home)
    moon.compute(home)

    target_timezone = pytz.timezone(selected_timezone)
    
    
    rising = home.previous_rising(sun).datetime().replace(tzinfo=pytz.utc).astimezone(target_timezone)
    sunrise = f'Sunrise is at {rising.hour}:{rising.minute}:{rising.second}'

    transit = home.next_transit(sun).datetime().replace(tzinfo=pytz.utc).astimezone(target_timezone)
    local_noon = f'Local noon is at {transit.hour}:{transit.minute}:{transit.second}'

    setting = home.next_setting(sun).datetime().replace(tzinfo=pytz.utc).astimezone(target_timezone)
    sunset = f'Sunset is at {setting.hour}:{setting.minute}:{setting.second}'

    moon_rising = home.previous_rising(moon).datetime().replace(tzinfo=pytz.utc).astimezone(target_timezone)
    moonrise = f'Moonrise is at {moon_rising.hour}:{moon_rising.minute}:{moon_rising.second}'

    moon_transit = home.next_transit(moon).datetime().replace(tzinfo=pytz.utc).astimezone(target_timezone)
    moon_transit_time = f'Moon transit is at {moon_transit.hour}:{moon_transit.minute}:{moon_transit.second}'

    moon_setting = home.next_setting(moon).datetime().replace(tzinfo=pytz.utc).astimezone(target_timezone)
    moonset = f'Moonset is at {moon_setting.hour}:{moon_setting.minute}:{moon_setting.second}'

    moon_phase = f'Moon phase is {moon.phase}%'

    day_length_sec = (home.next_setting(sun) - home.next_rising(sun)) * 86400
    hours = int(day_length_sec // 3600)
    minutes = int((day_length_sec % 3600) // 60)
    seconds = int(day_length_sec % 60)
    day_length = f"Day length is {hours} hours, {minutes} minutes, and {seconds} seconds"

    return sunrise, local_noon, sunset, moonrise, moon_transit_time, moonset, moon_phase, day_length

@app.route('/sunrise', methods=['GET', 'POST'])
def sunrise():
    time_zones = {'Asia/Kolkata': 'IST', 'America/New_York': 'EST', 'UTC': 'UTC'}

    date = dt.datetime.today().strftime('%Y-%m-%d') + ' 00:00:00'
    city = 'Ahmedabad'
    selected_timezone = 'Asia/Kolkata'

    if request.method == 'POST':
        date = request.form['date'] + ' 00:00:00'
        city = request.form['city']
        selected_timezone = request.form['timezone']

    result = calculate_astronomical_data(date, city, selected_timezone)
    if result is None:
        return render_template('error.html', message='Geocoding failed for the provided city.')

    sunrise, local_noon, sunset, moonrise, moon_transit_time, moonset, moon_phase, day_length = result

    return render_template('sunrise.html', date=date.split(' ')[0], city=city, selected_timezone=selected_timezone, sunrise=sunrise, local_noon=local_noon, sunset=sunset, moonrise=moonrise, moon_transit_time=moon_transit_time, moonset=moonset, moon_phase=moon_phase, day_length=day_length, time_zones=time_zones)




TIME_ZONES = ["Asia/Kolkata","UTC", "America/New_York", "Europe/London", "Asia/Tokyo", "Australia/Sydney", "Europe/Berlin"]  # Add more as needed

@app.route('/season')
def season():
    current_year = datetime.now().year
    return render_template('season.html', year=current_year, timezone="Asia/Kolkata", timezones=TIME_ZONES, results=calculate_for_year())

@app.route('/calculate', methods=['POST'])
def calculate():
    year = int(request.form['year'])
    timezone = request.form['timezone']
    return render_template('season.html', year=year, timezone=timezone, timezones=TIME_ZONES, results=calculate_for_year(year, timezone))

def calculate_for_year(year=None, timezone="UTC"):
    if year is None:
        year = datetime.now().year
    tz = pytz.timezone(timezone)

    ts = api.load.timescale()
    eph = api.load('de421.bsp')

    t0 = ts.utc(year, 1, 1)
    t1 = ts.utc(year, 12, 31)
    t, y = almanac.find_discrete(t0, t1, almanac.seasons(eph))

    results = []
    for yi, ti in zip(y, t):
        local_time = ti.astimezone(tz)
        formatted_time = local_time.strftime('%d %B %Y, %H:%M')
        event = f"{almanac.SEASON_EVENTS[yi]}: {formatted_time}"
        results.append(event)

    return results



@app.route("/weight", methods=["GET", "POST"])
def weight():
    planets = [
        "Earth", "Mercury", "Venus", "Sun", "Mars",
        "Jupiter", "Saturn", "Uranus", "Neptune"
    ]
    weights = {}

    if request.method == "POST":
        earth_weight = float(request.form.get("weight"))
        factors = [1, 0.38, 0.91, 27.9, 0.38, 2.34, 1.06, 0.92, 1.19] # Including Sun's gravitational factor
        weights = {planet: earth_weight * factor for planet, factor in zip(planets, factors)}

    return render_template('weight.html', weights=weights, planets=planets)






@app.route('/moon_phase', methods=['GET', 'POST'])
def moon_phase():
    date_str = None
    time_zone = "Asia/Kolkata"
    selected_timezones = ["Asia/Kolkata", "UTC", "America/New_York", "Europe/London", "Australia/Sydney"]

    if request.method == 'POST':
        date_str = request.form.get('date')
        time_zone = request.form.get('timezone')
        year, month, day = map(int, date_str.split('-'))
    else:
        now = datetime.utcnow()
        year, month, day = now.year, now.month, now.day
        date_str = now.strftime('%Y-%m-%d')

    tz = pytz.timezone(time_zone)  # Use pytz.timezone to get the timezone object
    selected_date = tz.localize(datetime(year, month, day))
    next_new_moon = ephem.localtime(ephem.next_new_moon(selected_date)).strftime('%Y-%b-%d')
    next_full_moon = ephem.localtime(ephem.next_full_moon(selected_date)).strftime('%Y-%b-%d')

    observer = ephem.Observer()
    observer.date = selected_date
    moon_phase_obj = ephem.Moon(observer)
    moon_phase = moon_phase_obj.phase
    percent_illumination = round(moon_phase, 2)

    new_moon_date_utc = ephem.previous_new_moon(selected_date)
    new_moon_date_local = ephem.localtime(new_moon_date_utc)
    new_moon_date_local_tz = tz.localize(new_moon_date_local)
    age_of_moon = round((selected_date - new_moon_date_local_tz).days + (selected_date - new_moon_date_local_tz).seconds / 86400, 2)

    if 0 <= moon_phase < 1.5:
        phase_name = 'New Moon'
    elif 1.5 <= moon_phase < 49.5:
        phase_name = 'Waxing Crescent'
    elif 49.5 <= moon_phase < 50.5:
        phase_name = 'First Quarter'
    elif 50.5 <= moon_phase < 99.5:
        phase_name = 'Waxing Gibbous'
    elif 99.5 <= moon_phase < 100.5:
        phase_name = 'Full Moon'
    elif 100.5 <= moon_phase < 149.5:
        phase_name = 'Waning Gibbous'
    elif 149.5 <= moon_phase < 150.5:
        phase_name = 'Last Quarter'
    else:
        phase_name = 'Waning Crescent'

    return render_template('moon_phase.html', percent=percent_illumination, phase_name=phase_name, date_str=date_str, next_new_moon_date=next_new_moon, next_full_moon_date=next_full_moon, timezones=selected_timezones, selected_timezone=time_zone, age_of_moon=age_of_moon)






ORBITAL_PERIODS = {
    'Earth': 365.25,
    'Mercury': 87.97,
    'Venus': 224.7,
    'Mars': 687,
    'Jupiter': 4332.59,
    'Saturn': 10759.22,
    'Uranus': 30688.5,
    'Neptune': 60182,
}

def calculate_age_on_planet(birth_date, orbital_period_in_earth_days):
    current_time = datetime.now()
    age_in_days = (current_time - birth_date).days
    age_on_planet = age_in_days / orbital_period_in_earth_days
    return round(age_on_planet, 1)  # Round to one decimal place

@app.route('/planetage', methods=['GET', 'POST'])
def planetage():
    results = {}
    birthdate = ""
    if request.method == 'POST':
        birthdate = request.form['birthdate']
        birth_date = datetime.strptime(birthdate, '%Y-%m-%d')
        for planet, period in ORBITAL_PERIODS.items():
            age = calculate_age_on_planet(birth_date, period)
            results[planet] = age

    return render_template('planetage.html', results=results, birthdate=birthdate)



@app.route('/julian', methods=['GET', 'POST'])
def julian():
    if request.method == 'POST':
        calendar_date = request.form['calendar_date']
        julian_date = convert_to_julian(calendar_date)
        return render_template('julian.html', calendar_date=calendar_date, julian_date=julian_date)
    return render_template('julian.html')

def convert_to_julian(calendar_date):
    try:
        year, month, day = map(int, calendar_date.split('-'))

        a = (14 - month) // 12
        y = year + 4800 - a
        m = month + 12 * a - 3

        julian_date = day + ((153 * m + 2) // 5) + (365 * y) + (y // 4) - (y // 100) + (y // 400) - 32045

        return str(julian_date)
    except Exception as e:
        return f"Error: {e}"
 
 
 
 
 
def generate_plot(city):
    location = geocoder.osm(city)

    if location.ok:
        latitude = location.latlng[0]
        longitude = location.latlng[1]
    else:
        return None, []

    stations_url = 'http://celestrak.com/NORAD/elements/stations.txt'
    tle_new_url = 'http://celestrak.com/NORAD/elements/tle-new.txt'
    satellites = load.tle_file(stations_url)
    satellites.extend(load.tle_file(tle_new_url))

    ts = load.timescale()
    t = ts.now()

    observer = Topos(latitude_degrees=latitude, longitude_degrees=longitude)

    altitudes = []
    azimuths = []
    names = []

    for satellite in satellites:
        difference = satellite - observer
        topocentric = difference.at(t)
        
        alt, az, distance = topocentric.altaz()
        
        if alt.degrees > 0:
            altitudes.append(90 - alt.degrees)
            azimuths.append(az.degrees)
            names.append(satellite.name)

    # Create a subplot grid for better layout
    fig, ax = plt.subplots(figsize=(10, 15), subplot_kw={'projection': 'polar'})
    ax.set_facecolor('black')
    
    # Plot satellites with markers, arrows, and labels
    for i, (azimuth, altitude, name) in enumerate(zip(azimuths, altitudes, names)):
        color = plt.cm.rainbow(i / len(azimuths))
        marker = ax.plot(np.radians(azimuth), altitude, marker='o', markersize=8, color=color, label=name)[0]
        
        # Add an arrow to indicate the satellite's direction of motion
        arrow_length = 0.1
        arrow_dx = arrow_length * np.sin(np.radians(azimuth))
        arrow_dy = arrow_length * np.cos(np.radians(azimuth))
        ax.arrow(np.radians(azimuth), altitude, arrow_dx, arrow_dy, head_width=0.05, head_length=0.1, fc=color, ec=color)
        
        # Set a label for the marker
        marker.set_label(name)
    
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_yticklabels([])
    plt.xticks(fontsize=20)
    plt.title("Visible Satellites", fontsize=20, fontweight='bold')
    
    # Adjust legend placement
    legend = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2)
    
    # Add cardinal direction labels at the bottom
    directions = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
    angles = [0, 45, 90, 135, 180, 225, 270, 315]
    ax.set_xticks(np.radians(angles))
    ax.set_xticklabels(directions)
    
    buffer = io.BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)
    image_base64 = base64.b64encode(buffer.read()).decode()
    buffer.close()

    satellite_details = []
    for name, alt, az in zip(names, altitudes, azimuths):
        satellite_details.append({'name': name, 'altitude': round(alt, 3), 'azimuth': round(az, 3)})

    return image_base64, satellite_details

# Generate satellite details and plot for Ahmedabad when the app starts
default_city = 'Ahmedabad'  # Replace with the desired default city
default_location = geocoder.osm(default_city)

if default_location.ok:
    default_latitude = default_location.latlng[0]
    default_longitude = default_location.latlng[1]
    default_image, default_satellite_details = generate_plot(default_city)
else:
    default_image = None
    default_satellite_details = []

@app.route('/satellite_track', methods=['GET', 'POST'])
def satellite_track():
    # Initialize variables to None
    image = None
    satellite_details = []
    entered_city = ''

    if request.method == 'POST':
        # Get the entered city value
        entered_city = request.form.get('city', '')

        if entered_city:
            # Generate plot for the selected city
            image, satellite_details = generate_plot(entered_city)
    else:
        # If no city is entered in the form, generate plot for a default city
        entered_city = 'Ahmedabad'  # Default city
        image, satellite_details = generate_plot(entered_city)

    return render_template('satellite_track.html', image=image, satellite_details=satellite_details, entered_city=entered_city)
    
    
    
def read_csv(filename):
    data = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        headers = lines[0].strip().split(',')
        for line in lines[1:]:
            values = line.strip().split(',')
            entry = dict(zip(headers, values))
            data.append(entry)
    return data

solar_eclipse_data = read_csv('Solar_e.csv')
lunar_eclipse_data = read_csv('lunar_e.csv')

@app.route('/lunar', methods=['GET', 'POST'])
def lunar():
    if request.method == 'POST':
        year = int(request.form['year'])
    else:
        year = 2023  # Default year
    
    solar_eclipses = [entry for entry in solar_eclipse_data if entry['Year'] == str(year)]
    lunar_eclipses = [entry for entry in lunar_eclipse_data if entry['Year'] == str(year)]
    
    return render_template('lunar.html', solar_eclipses=solar_eclipses, lunar_eclipses=lunar_eclipses, year=year)




@app.route('/solar_system')
def index():
    return render_template('solar_system.html')

@app.route('/get_orbits', methods=['POST'])
def get_orbits():
    request_data = request.get_json()
    date = request_data['date']
    time = request_data['time']
    t = Time(date + "T" + time, scale="tdb")
    planets = ['mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']
    data = []

    for planet in planets:
        pos, vel = get_body_barycentric_posvel(planet, t)
        
        # Calculate speed in AU per day
        speed_au_per_day = np.sqrt(vel.x.value**2 + vel.y.value**2 + vel.z.value**2)
        
        # Convert speed to m/s
        speed_m_per_s = speed_au_per_day * 149597870.7 * 1000 / (24 * 3600)
        
        theta = np.linspace(0, 2 * np.pi, 1000)
        r = np.sqrt(pos.x.value**2 + pos.y.value**2)
        x_orbit = (r * np.cos(theta)).tolist()
        y_orbit = (r * np.sin(theta)).tolist()

        data.append({
            'planet': planet.capitalize(),
            'x_orbit': x_orbit,
            'y_orbit': y_orbit,
            'x_pos': pos.x.value,
            'y_pos': pos.y.value,
            'distance_from_sun': r,
            'speed': speed_m_per_s  # Speed in m/s
        })

    return jsonify(data)
 

timezones = pytz.all_timezones


@app.route('/time_zone', methods=['GET', 'POST'])
def time_zone():
    input_time = ''
    from_timezone = ''
    to_timezone = ''
    converted_datetime = ''
    calculated_timezone = ''

    if request.method == 'POST':
        input_time = request.form['input_time']
        from_timezone = request.form['from_timezone']
        to_timezone = request.form['to_timezone']

        input_datetime = datetime.strptime(input_time, '%Y-%m-%dT%H:%M')

        from_zone = pytz.timezone(from_timezone)
        input_datetime = from_zone.localize(input_datetime)

        to_zone = pytz.timezone(to_timezone)
        converted_datetime = input_datetime.astimezone(to_zone)
        calculated_timezone = to_timezone

    return render_template('time_zone.html', timezones=timezones, input_time=input_time,
                           from_timezone=from_timezone, to_timezone=to_timezone,
                           converted_datetime=converted_datetime,
                           calculated_timezone=calculated_timezone)



def dms_to_dd(degrees, minutes, seconds):
    dd = degrees + minutes / 60 + seconds / 3600
    return dd

@app.route('/dsmdd', methods=['GET', 'POST'])
def dsmdd():
    latitude = longitude = None
    lat_degrees = lat_minutes = lat_seconds = None
    lon_degrees = lon_minutes = lon_seconds = None
    
    if request.method == 'POST':
        lat_degrees = float(request.form['lat_degrees'])
        lat_minutes = float(request.form['lat_minutes'])
        lat_seconds = float(request.form['lat_seconds'])

        lon_degrees = float(request.form['lon_degrees'])
        lon_minutes = float(request.form['lon_minutes'])
        lon_seconds = float(request.form['lon_seconds'])

        latitude = dms_to_dd(lat_degrees, lat_minutes, lat_seconds)
        longitude = dms_to_dd(lon_degrees, lon_minutes, lon_seconds)

    return render_template('dsmdd.html', latitude=latitude, longitude=longitude,
                           lat_degrees=lat_degrees, lat_minutes=lat_minutes, lat_seconds=lat_seconds,
                           lon_degrees=lon_degrees, lon_minutes=lon_minutes, lon_seconds=lon_seconds)




units = {
    'light_year': 'Light Years',
    'parsec': 'Parsecs',
    'kilometer': 'Kilometers',
    'mile': 'Miles',
    'au': 'Astronomical Units'
}

@app.route('/lightyear', methods=['GET', 'POST'])
def lightyear():
    from_unit = ''
    to_unit = ''
    value = ''
    converted_value = ''

    if request.method == 'POST':
        from_unit = request.form['from_unit']
        to_unit = request.form['to_unit']
        value = float(request.form['value'])

        if from_unit != to_unit:
            converted_value = convert_distance(value, from_unit, to_unit)

    return render_template('lightyear.html', units=units, from_unit=from_unit, to_unit=to_unit,
                           value=value, converted_value=converted_value)

def convert_distance(value, from_unit, to_unit):
    conversions = {
        ('light_year', 'parsec'): lambda x: x * 3.262,
        ('light_year', 'kilometer'): lambda x: x * 9.461e+12,
        ('light_year', 'mile'): lambda x: x * 5.8786254e+12,
        ('light_year', 'au'): lambda x: x * 63241.077084,
        ('parsec', 'light_year'): lambda x: x / 3.262,
        ('parsec', 'kilometer'): lambda x: x * 3.086e+13,
        ('parsec', 'mile'): lambda x: x * 1.9173135e+13,
        ('parsec', 'au'): lambda x: x * 206265,
        ('kilometer', 'light_year'): lambda x: x / 9.461e+12,
        ('kilometer', 'parsec'): lambda x: x / 3.086e+13,
        ('kilometer', 'mile'): lambda x: x * 0.621371,
        ('kilometer', 'au'): lambda x: x / 149597870.7,
        ('mile', 'light_year'): lambda x: x / 5.8786254e+12,
        ('mile', 'parsec'): lambda x: x / 1.9173135e+13,
        ('mile', 'kilometer'): lambda x: x / 0.621371,
        ('mile', 'au'): lambda x: x / 92955807.3,
        ('au', 'light_year'): lambda x: x / 63241.077084,
        ('au', 'parsec'): lambda x: x / 206265,
        ('au', 'kilometer'): lambda x: x * 149597870.7,
        ('au', 'mile'): lambda x: x * 92955807.3
    }

    return conversions[(from_unit, to_unit)](value)


### Twilight



important_timezones = ['Asia/Kolkata', 'UTC', 'US/Pacific', 'Europe/London', 'Asia/Shanghai', 'Australia/Sydney']

def get_twilight_times(city_name, date_input, tz):
    location = geocoder.osm(city_name)
    lat, lon = location.latlng

    date_input = date_input or dt.datetime.now().date().isoformat()
    year, month, day = map(int, date_input.split('-'))
    date = dt.datetime(year, month, day)

    zone = timezone(tz)
    date = zone.localize(date)

    midnight = date.replace(hour=0, minute=0, second=0, microsecond=0)
    next_midnight = midnight + dt.timedelta(days=1)

    ts = load.timescale()
    t0 = ts.from_datetime(midnight)
    t1 = ts.from_datetime(next_midnight)
    eph = load('de421.bsp')

    observer = wgs84.latlon(lat * N, lon * E)

    f = almanac.dark_twilight_day(eph, observer)
    times, events = almanac.find_discrete(t0, t1, f)

    previous_e = f(t0).item()
    results = []
    for t, e in zip(times, events):
        tstr = t.astimezone(zone).time().strftime('%I:%M %p')
        if previous_e < e:
            results.append(f'{almanac.TWILIGHTS[e]} starts at: {tstr}')
        else:
            results.append(f'{almanac.TWILIGHTS[previous_e]} ends at: {tstr}')
        previous_e = e

    return results

@app.route('/twilight', methods=['GET', 'POST'])
def twilight():
    city_name = request.form.get('city', 'Mumbai')
    date_input = request.form.get('date', dt.datetime.now().date().isoformat())
    tz = request.form.get('timezone', 'Asia/Kolkata')

    results = get_twilight_times(city_name, date_input, tz)

    return render_template('twilight.html', results=results, timezones=important_timezones, date_input=date_input)


## GPS Time 

def convert_to_gps_time(calendar_time):
    gps_epoch = datetime(1980, 1, 6, 0, 0, 0)
    time_difference = calendar_time - gps_epoch
    gps_seconds = int(time_difference.total_seconds())
    return gps_seconds

@app.route('/gps_time', methods=['GET', 'POST'])
def gps_time():
    gps_time = None

    if request.method == 'POST':
        input_datetime = request.form['input_datetime']
        try:
            calendar_time = datetime.strptime(input_datetime, '%Y-%m-%dT%H:%M')
            gps_time = convert_to_gps_time(calendar_time)
        except ValueError:
            error_message = "Invalid datetime format. Please use YYYY-MM-DDTHH:MM."
            return render_template('index.html', error_message=error_message)

    return render_template('gps_time.html', gps_time=gps_time)



## Solar Eclipse



@app.route('/solar_eclipse', methods=['GET', 'POST'])
def solar_eclipse():
    results = []
    year_entered = ""  # Default value is an empty string

    if request.method == 'POST':
        year_entered = request.form['year']
        results = get_eclipses(int(year_entered))
    
    return render_template('solar_eclipse.html', results=results, year_entered=year_entered)

def get_eclipses(year):
    curtime = datetime(year, 1, 1, 0, 0, 0)
    endtime = datetime(year, 12, 31, 23, 59, 59)
    moon = ephem.Moon()
    sun = ephem.Sun()
    observer = ephem.Observer()
    observer.elevation = -6371000
    observer.pressure = 0
    
    eclipse_dates = []

    while curtime <= endtime:
        observer.date = curtime.strftime('%Y/%m/%d %H:%M:%S')
        moon.compute(observer)
        sun.compute(observer)
        sep = abs((float(ephem.separation(moon, sun)) / 0.01745329252))
        if sep < 1.59754941:
            eclipse_dates.append(curtime.strftime('%B %d, %Y'))
            curtime += timedelta(days=1)
        else:
            curtime += timedelta(minutes=5)

    return eclipse_dates


#######################  CONJUNCTION



@app.route("/planet_conjunction", methods=["GET", "POST"])
def planet_conjunction():
    year = ""  # Initialize the year variable
    if request.method == "POST":
        year = request.form.get("year")  # Get the entered year
        target_timezone = request.form.get("timezone")

        data = []
        headings = None

        csv_dir = "csv_files2"
        for csv_file in os.listdir(csv_dir):
            if csv_file.endswith(".csv"):
                with open(os.path.join(csv_dir, csv_file), "r") as csvfile:
                    reader = csv.reader(csvfile)
                    for row in reader:
                        if headings is None:
                            headings = row
                        if len(row) > 0 and row[1].startswith(year):
                            data.append(row)

        return render_template("planet_conjunction.html", headings=headings, data=data, target_timezone=target_timezone, year=year)

    return render_template("planet_conjunction.html", headings=None, data=None, year=year)  # Pass the year to the template

####################


@app.route("/moon_conj", methods=["GET", "POST"])
def moon_conj():
    year = ""
    if request.method == "POST":
        year = request.form.get("year")
        target_timezone = request.form.get("timezone")

        data = []
        headings = None

        csv_dir = "csv_files3"
        for csv_file in os.listdir(csv_dir):
            if csv_file.endswith(".csv"):
                with open(os.path.join(csv_dir, csv_file), "r") as csvfile:
                    reader = csv.reader(csvfile)
                    for row in reader:
                        if headings is None:
                            headings = row
                        if len(row) > 0 and row[1].startswith(year):
                            data.append(row)

        return render_template("moon_conj.html", headings=headings, data=data, target_timezone=target_timezone, year=year)

    return render_template("moon_conj.html", headings=None, data=None, year=year)




def haversine_distance(lat1, lon1, lat2, lon2):
    R = 6371.0
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)
    d_lat = lat2_rad - lat1_rad
    d_lon = lon2_rad - lon1_rad
    a = math.sin(d_lat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(d_lon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = R * c
    return distance

@app.route('/latlon_distance', methods=['GET', 'POST'])
def latlon_distance():
    lat1 = lat2 = lon1 = lon2 = distance = None
    if request.method == 'POST':
        lat1 = float(request.form['lat1'])
        lon1 = float(request.form['lon1'])
        lat2 = float(request.form['lat2'])
        lon2 = float(request.form['lon2'])
        distance = haversine_distance(lat1, lon1, lat2, lon2)
    return render_template('latlon_distance.html', lat1=lat1, lon1=lon1, lat2=lat2, lon2=lon2, distance=distance)







#################################### map folium dark place#####################

dark_sky_locations = [
    {"name": "Spiti Valley", "latitude": 32.2708, "longitude": 78.0618},
    {"name": "Nubra Valley / Ladakh", "latitude": 34.1526, "longitude": 77.5770},
    {"name": "Kutch Desert, Gujarat", "latitude": 23.7337, "longitude": 69.8597},
    {"name": "B R Hills, Karnataka", "latitude": 11.9287, "longitude": 77.1195},
    {"name": "Kanha National Park, MP", "latitude": 22.3309, "longitude": 80.6110},
    {"name": "Tal Chhapar Sanctuary, Rajasthan", "latitude": 27.0074, "longitude": 74.9995},
    {"name": "Valley of Flowers National Park, Uttarakhand", "latitude": 30.7272, "longitude": 79.6062},
    {"name": "Agumbe, Karnataka", "latitude": 13.5073, "longitude": 75.0831},
    {"name": "Sariska Tiger Reserve, Rajasthan", "latitude": 27.1975, "longitude": 76.3669},
    {"name": "Matheran, Maharashtra", "latitude": 18.9865, "longitude": 73.2659},
    {"name": "Mount Abu, Rajasthan", "latitude": 24.5925, "longitude": 72.7156},
    {"name": "Chandratal Lake, Himachal Pradesh", "latitude": 32.4545, "longitude": 77.6152},
    {"name": "Sundarbans National Park, West Bengal", "latitude": 21.9497, "longitude": 88.9002},
    {"name": "Gir Forest National Park, Gujarat", "latitude": 21.1897, "longitude": 70.9421},
    {"name": "Ganeshgudi, Karnataka", "latitude": 15.1561, "longitude": 74.5007},
    {"name": "Dhanaulti, Uttarakhand", "latitude": 30.4221, "longitude": 78.2935},
    {"name": "Kausani, Uttarakhand", "latitude": 29.8456, "longitude": 79.6012},
    {"name": "Vattakanal, Tamil Nadu", "latitude": 10.2397, "longitude": 77.4844},
    {"name": "Pelling, Sikkim", "latitude": 27.3170, "longitude": 88.6254},
    {"name": "Araku Valley, Andhra Pradesh", "latitude": 18.3344, "longitude": 82.8687},
    {"name": "Cherrapunji, Meghalaya", "latitude": 25.2966, "longitude": 91.5822},
    {"name": "Pachmarhi, Madhya Pradesh", "latitude": 22.4674, "longitude": 78.4342},
    {"name": "Binsar, Uttarakhand", "latitude": 29.6621, "longitude": 79.6555},
    {"name": "Ramgarh, Uttarakhand", "latitude": 29.5134, "longitude": 79.3241},
    # Add more locations here
    # ...
]

@app.route('/darkmap')
def darkmap():
    m = folium.Map(location=[20.5937, 78.9629], zoom_start=5.4)

    for location in dark_sky_locations:
        folium.Marker(
            [location["latitude"], location["longitude"]],
            popup=location["name"],
            icon=folium.Icon(color='black')
        ).add_to(m)

    return m._repr_html_()


@app.route('/spacemap')
def spacemap():
    # Create a map centered at a specific location
    map_center = [20, 0]  # Adjust the center as needed
    mymap = folium.Map(location=map_center, zoom_start=3)

    agencies = [
        ("NASA", [38.8838, -77.0164]),
        ("Roscosmos", [55.7058, 37.6466]),
        ("ESA", [48.8534, 2.3488]),
        ("CNSA", [40.0799, 116.6031]),
        ("ISRO", [12.9716, 77.5946]),
        ("JAXA", [35.6895, 139.6917]),
        ("CSA", [45.5007, -73.4162]),
        ("UK Space Agency", [51.5685, -1.7722]),
        ("UAE Space Agency", [24.4539, 54.3773]),
        ("DLR", [50.9350, 6.9385]),
        ("ASI", [41.9028, 12.4964]),
        ("KARI", [36.3504, 127.3848]),
        ("CNES", [48.8415, 2.3036]),
        ("CONAE", [-31.4201, -64.1888]),
        ("SUPARCO", [33.6844, 73.0479]),
        ("ANGKASA", [3.1390, 101.6869]),
        ("ISA", [32.0853, 34.7818]),
        ("ROSA", [44.4268, 26.1025]),
        ("CAST", [40.0439, 116.3215]),
        ("CNPE", [-23.2237, -45.9009]),
        ("NSAU", [50.4430, 30.5220]),
        ("ASIIN", [50.9020, 7.0955]),
        ("COSMO-SkyMed", [41.8106, 12.5213]),
        ("ROSCOMOS", [54.9011, 83.0000]),
        # Add more agencies and coordinates here...
    ]

   
    
    for agency, coords in agencies:
        folium.Marker(location=coords, popup=agency).add_to(mymap)

    # Save the map as an HTML string
    map_html = mymap._repr_html_()

    return render_template('spacemap.html', map_html=map_html)




launch_pads = [
    {"name": "Kennedy Space Center LC-39A", "lat": 28.5721, "lon": -80.6480},
    {"name": "Baikonur Cosmodrome Gagarin's Start", "lat": 45.9654, "lon": 63.3058},
    {"name": "Guiana Space Centre Ariane Launch Complex", "lat": 5.2360, "lon": -52.7763},
    {"name": "Satish Dhawan Space Centre Second Launch Pad", "lat": 13.7337, "lon": 80.2351},
    {"name": "Jiuquan Satellite Launch Center LC-43", "lat": 40.9585, "lon": 100.2917},
    {"name": "Vandenberg SLC-4E", "lat": 34.6400, "lon": -120.6100},
    {"name": "Tanegashima Yoshinobu Launch Complex", "lat": 30.3850, "lon": 130.9656},
    {"name": "Plesetsk LC-43/4", "lat": 62.9271, "lon": 40.5785},
    {"name": "Wenchang Spacecraft Launch Site", "lat": 19.6145, "lon": 110.9513},
    {"name": "Yasny Cosmodrome", "lat": 51.0250, "lon": 59.8600},
    {"name": "Naro Space Center", "lat": 34.3903, "lon": 127.3654},
    {"name": "Andøya Space Center", "lat": 69.2924, "lon": 16.0299},
    {"name": "Kodiak Launch Complex", "lat": 57.4700, "lon": -152.3080},
    {"name": "Dombarovsky Air Base Yasny Cosmodrome", "lat": 51.0250, "lon": 59.8600},
    {"name": "Alcântara Spaceport CLA-41", "lat": -2.3750, "lon": -44.3875},
    {"name": "Mid-Atlantic Regional Spaceport LP-0B", "lat": 37.8330, "lon": -75.4880},
    {"name": "Semnan Spaceport", "lat": 35.2344, "lon": 53.9369},
    {"name": "Uchinoura Space Center", "lat": 31.2510, "lon": 131.0843},
    {"name": "Vostochny Cosmodrome", "lat": 51.8844, "lon": 128.3331},
    {"name": "Jiuquan LC-4", "lat": 40.9585, "lon": 100.2917},
    {"name": "Spaceport America", "lat": 32.9903, "lon": -106.9744},
    {"name": "Yasny Cosmodrome", "lat": 51.0250, "lon": 59.8600},
    {"name": "North Korea's Sohae Satellite Launching Station", "lat": 39.6603, "lon": 124.7059},
    # Add more launch pads with their latitudes and longitudes
    # ...
]

@app.route('/launchpads')
def launchpads():
    launch_map = folium.Map(location=[20, 0], zoom_start=3)
    
    for pad in launch_pads:
        folium.Marker(
            [pad["lat"], pad["lon"]],
            popup=pad["name"],
            icon=folium.Icon(icon="rocket", prefix="fa", color="red")
        ).add_to(launch_map)
    
    return render_template('launchpads.html', launch_map=launch_map._repr_html_())



def degrees_to_hms(degrees):
    total_seconds = int(degrees * 3600)
    hours = total_seconds // 3600
    remaining_seconds = total_seconds % 3600
    minutes = remaining_seconds // 60
    seconds = remaining_seconds % 60
    return hours, minutes, seconds



@app.route('/degree_convert', methods=['GET', 'POST'])
def degree_convert():
    result = None
    degrees_value = ""
    explanation = (
        "In astronomy, celestial coordinates often use degrees, but for visualizing the night sky, "
        "hours, minutes, and seconds are more intuitive. This converter helps translate between them."
    )
    if request.method == 'POST':
        degrees = float(request.form['degrees'])
        hours, minutes, seconds = degrees_to_hms(degrees)
        result = f"{degrees} degrees is equivalent to {hours} hours, {minutes} minutes, {seconds} seconds"
        degrees_value = degrees
    return render_template('degree_convert.html', result=result, explanation=explanation, degrees_value=degrees_value)





############ STAR



from datetime import datetime
from pytz import timezone
from starplot import ZenithPlot
from starplot.styles import PlotStyle, extensions
import geocoder
import io
import base64
import matplotlib


@app.route("/starchart", methods=["GET", "POST"])
def starchart():
    # Default values
    default_city_name = "Ahmedabad"
    default_timezone = "Asia/Kolkata"

    # Get the current date and time
    current_datetime = datetime.now().astimezone(timezone(default_timezone))
    default_date = current_datetime.strftime("%Y-%m-%d")
    default_time = current_datetime.strftime("%H:%M")

    if request.method == "POST":
        city_name = request.form.get("city_name")
        date_str = request.form.get("date")
        time_str = request.form.get("time")
        timezone_name = request.form.get("timezone")

        try:
            location = geocode_city(city_name)

            lat = location.latlng[0]
            lon = location.latlng[1]

            tz = timezone(timezone_name)
            dt_str = f"{date_str} {time_str}"
            dt = datetime.strptime(dt_str, "%Y-%m-%d %H:%M")
            dt = tz.localize(dt)

            p = ZenithPlot(
                lat=lat,
                lon=lon,
                dt=dt,
                limiting_magnitude=4.6,
                style=PlotStyle().extend(
                    extensions.BLUE_MEDIUM,
                    extensions.ZENITH,
                ),
                resolution=2000,
                adjust_text=True,
            )

            # Create a BytesIO object to store the image data
            img_buffer = io.BytesIO()
            p.export(img_buffer, format="png")
            img_buffer.seek(0)  # Reset buffer position to the beginning

            # Encode the image data as Base64
            img_data_base64 = base64.b64encode(img_buffer.read()).decode("utf-8")

            return render_template(
                "starchart.html",
                image_data=img_data_base64,
                message=f"Star chart for {city_name} on {date_str} at: {time_str}",
                city_name=city_name,
                date=date_str,
                time=time_str,
                timezone=timezone_name,
            )
        except Exception as e:
            return render_template(
                "starchart.html",
                image_data=None,
                message=f"Error: {str(e)}",
            )

    # Handle the case when no specific city is selected
    return render_template(
        "starchart.html",
        image_data=None,
        message="Please select a city and provide date and time information.",
        city_name="",
        date=default_date,
        time=default_time,
        timezone=default_timezone,
    )

def geocode_city(city_name):
    try:
        location = geocoder.osm(city_name)
        return location
    except Exception as e:
        raise Exception("Failed to retrieve location data")
    





###############Sat Orbit

@app.route('/sat_orbit', methods=['GET', 'POST'])
def sat_orbit():
    if request.method == 'POST':
        try:
            altitude_km = float(request.form['altitude'])
            altitude_m = altitude_km * 1000  # Convert kilometers to meters
            radius = altitude_m + 6371000  # Earth's radius is 6371 km in meters
            gravitational_constant = 6.67430e-11
            earth_mass = 5.972e24

            orbital_period = 2 * math.pi * math.sqrt((radius ** 3) / (gravitational_constant * earth_mass))
            
            # Format the result with up to two decimal places
            formatted_result = "{:.2f}".format(orbital_period / 3600)  # Convert seconds to hours
            
            return render_template('sat_orbit.html', result=formatted_result)
        except ValueError:
            error_message = "Invalid input. Please enter a numeric value for altitude in kilometers."
            return render_template('sat_orbit.html', error=error_message)

    return render_template('sat_orbit.html')



########### Apogee


@app.route('/apogee', methods=['GET', 'POST'])
def apogee():
    eccentricity, semimajor_axis, ha, hp = None, None, None, None
    
    if request.method == 'POST':
        eccentricity = request.form['eccentricity']
        semimajor_axis = request.form['semimajor_axis']
        
        if eccentricity and semimajor_axis:
            eccentricity = float(eccentricity)
            semimajor_axis = float(semimajor_axis)
            
            # Calculation
            R = 6371  # Mean Earth's radius (km)
            ra = semimajor_axis * (1 + eccentricity)  # Radius Vector at apogee (km)
            rp = semimajor_axis * (1 - eccentricity)  # Radius Vector at perigee (km)
            ha = round(ra - R, 2)  # Apogee height (km)
            hp = round(rp - R, 2)  # Perigee height (km)


            hp = round(ra - R, 2)
            ra = semimajor_axis 
    


    return render_template('apogee.html', eccentricity=eccentricity, semimajor_axis=semimajor_axis, ha=ha, hp=hp)




# Define constants
G = 6.67430e-11  # Gravitational constant (m^3/kg/s^2)
M = 5.972e24     # Mass of the Earth (kg)

@app.route("/sat_speed", methods=["GET", "POST"])
def sat_speed():
    altitude = None
    speed = None

    if request.method == "POST":
        altitude = float(request.form["altitude"])  # Altitude in kilometers
        altitude_meters = altitude * 1000  # Convert altitude to meters

        # Calculate satellite speed in m/s using the formula
        speed = (G * M / (altitude_meters + 6371000)) ** 0.5

    return render_template("sat_speed.html", altitude=altitude, speed=speed)







astronauts_url = "http://api.open-notify.org/astros.json"
iss_location_url = "http://api.open-notify.org/iss-now.json"

# Function to fetch astronaut data and ISS location
def get_astronauts_and_iss_location():
    try:
        # Get astronaut data
        response_astronauts = requests.get(astronauts_url)
        astronauts_data = response_astronauts.json()
        
        # Get ISS location data
        response_iss_location = requests.get(iss_location_url)
        iss_location_data = response_iss_location.json()
        
        # Extract relevant information
        num_astronauts = astronauts_data["number"]
        astronauts = astronauts_data["people"]
        iss_lat = iss_location_data["iss_position"]["latitude"]
        iss_lon = iss_location_data["iss_position"]["longitude"]
        
        return num_astronauts, astronauts, iss_lat, iss_lon
        
    except requests.exceptions.RequestException as e:
        return None

# Define a route for the home page
@app.route('/iss')
def iss():
    astronaut_data = get_astronauts_and_iss_location()
    return render_template('iss.html', astronaut_data=astronaut_data)



if __name__ == '__main__':
    # Run the app in debug mode on port 5000
    app.run(debug=True, port=5000, host='0.0.0.0')
