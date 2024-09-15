f''' helpers.py
    Copyright (C) 2024 github.com/brandongrahamcobb

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

from bot.main import Lucy
from datetime import datetime
from discord.ext import commands
from gradio_client import Client
from io import BytesIO
from matplotlib import pyplot as plt
from os import makedirs
from os.path import abspath, dirname, exists, expanduser, isfile, join
from PIL import Image, ImageFont, ImageDraw
from random import randint
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, DataStructs, Draw, rdDepictor, rdFingerprintGenerator, rdFMCS
rdDepictor.SetPreferCoordGen(True)
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.DataStructs import FingerprintSimilarity, TanimotoSimilarity
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from typing import List
from webdriver_manager.chrome import ChromeDriverManager
from typing import List, Optional, Dict, Any

global logger

import asyncio
import colorsys
import datetime as dt
import discord
import emoji as emoji_lib
import itertools
import logging
import logging.handlers
import math
import openai
import os
import pubchempy as pcp
import random
import requests
import traceback
import unicodedata
import yaml

current_date = dt.datetime.now().strftime('%d%m%y')

home = expanduser('~')

path_base = join(home, 'Documents', 'src', 'lucy')
path_ai_cog = join(path_base, 'cogs', 'ai_cog.py')
path_config_yaml = join(home, '.config', 'lucy', 'config.yaml')
path_log = join(home, '.log', 'discord.log')
path_users_yaml = join(home, '.config', 'lucy', 'users.yaml')

def add_watermark(image: BytesIO, watermark_text: str = '~LUCY~') -> BytesIO:
    RGB_image = Image.open(image)
    RGBA_image = RGB_image.convert('RGBA')
    draw = ImageDraw.Draw(RGBA_image)
    width, height = RGBA_image.size
    diagonal = math.sqrt(width**2 + height**2)
    font_size = int(diagonal / 15)
    try:
        font = ImageFont.truetype('resources/Roboto-Regular.ttf', font_size)  # Replace with the path to your font file
    except IOError:
        font = ImageFont.load_default()
    while True:
        bbox = draw.textbbox((0, 0), watermark_text, font=font)
        text_width = bbox[2] - bbox[0]
        if text_width <= 512 or font_size <= 1:
            break
        font_size -= 1
        font = ImageFont.truetype('resources/Roboto-Regular.ttf', font_size)  # Replace with the path to your font file
    text_height = bbox[3] - bbox[1]
    text_x = (width - text_width) / 2
    text_y = height - (2 * text_height)
    watermark_image = Image.new('RGBA', RGBA_image.size, (0, 0, 0, 0))
    watermark_draw = ImageDraw.Draw(watermark_image)
    watermark_draw.text((text_x, text_y), watermark_text, font=font, fill=(255, 255, 255, 64))
    mask = watermark_image.split()[3]
    RGBA_image.paste(watermark_image, (0, 0), mask)
    output = BytesIO()
    RGBA_image.save(output, format= 'PNG')
    output.seek(0)
    return output

def adjust_hue_and_saturation(image, hue_shift, saturation_shift) -> BytesIO:
    image = image.convert('RGB')
    pixels = list(image.getdata())
    adjusted_pixels = []
    for r, g, b in pixels:
        h, s, v = colorsys.rgb_to_hsv(r / 255.0, g / 255.0, b / 255.0)
        h = (h + hue_shift / 360.0) % 1.0
        s = min(max(s + saturation_shift / 100.0, 0), 1)
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        adjusted_pixels.append((int(r * 255), int(g * 255), int(b * 255)))
    new_image = Image.new('RGB', image.size)
    new_image.putdata(adjusted_pixels)
    output = BytesIO()
    new_image.save(output, format= 'PNG')
    output.seek(0)
    return output

def combine(images: list, names: list) -> BytesIO:
    combined_images = []
    for index, (bytes_io, name) in enumerate(zip(images, names)):
        img = Image.open(bytes_io)
        inverted_image = Image.eval(img, lambda x: 255 - x)
        img_bytes = BytesIO()
        inverted_image.save(img_bytes, format='PNG')
        img_bytes_final = add_watermark(img_bytes, watermark_text=name)
        img_final = Image.open(img_bytes_final)
        combined_images.append(img_final)
    widths, heights = zip(*(img.size for img in combined_images))
    total_width = sum(widths)
    max_height = max(heights)
    combined_img = Image.new('RGB', (total_width, max_height))
    x_offset = 0
    for img in combined_images:
        combined_img.paste(img, (x_offset, 0))
        x_offset += img.width
    output = adjust_hue_and_saturation(combined_img, hue_shift=-180, saturation_shift=160)
    output.seek(0)
    return output

def draw_fingerprint(pair) -> BytesIO:
#    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, countSimulation=True)
#    def get_fp(mol, *args, **kwargs):
#        return mfpgen.GetFingerprint(mol)
    d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
    d2d.prepareMolsBeforeDrawing = False
    Options = d2d.drawOptions()
    Options.prepareMolsBeforeDrawing = False
    Options.includeMetadata = False
    Options.bondLineWidth = 4.0
    d2d.SetDrawOptions(Options)
    mol1 = rdMolDraw2D.PrepareMolForDrawing(pair[0], kekulize=True)
    mol1.UpdatePropertyCache(False)
    mol2 = rdMolDraw2D.PrepareMolForDrawing(pair[1], kekulize=True)
    mol2.UpdatePropertyCache(False)
    fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol1, mol2, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d, drawingOptions=Options)
#    fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(pair[1], pair[0], get_fp, draw2d=d2d) #colorMap=brighter_color, draw2d=d2d)
    d2d.FinishDrawing()
    drawing = d2d.GetDrawingText()
    output = BytesIO(drawing)
    return output

def draw_watermarked_molecule(molecule) -> BytesIO:
    resolved_name = get_molecule_name(molecule)
    d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
    Options = d2d.drawOptions()
    Options.prepareMolsBeforeDrawing = False
    Options.includeMetadata = False
    Options.bondLineWidth = 4.0
    d2d.SetDrawOptions(Options)
    rdMolDraw2D.SetDarkMode(Options)
    mol = rdMolDraw2D.PrepareMolForDrawing(molecule, kekulize=True)
    mol.UpdatePropertyCache(False)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    drawing = d2d.GetDrawingText()
    output = add_watermark(BytesIO(drawing), watermark_text=resolved_name)
    return output

def get_emoji(emoji_character):
    if emoji_character is None:
        return 'Please provide a Unicode emoji character.'
    unicode_name = emoji_lib.demojize(emoji_character)
    if unicode_name.startswith(':') and unicode_name.endswith(':'):
        unicode_name = unicode_name[1:-1]
        code_points = ' '.join(f'U+{ord(c):04X}' for c in emoji_character)
        description = unicodedata.name(emoji_character, 'No description available')
        unicode_info = (
            f'**Unicode Emoji Character:** {emoji_character}\n'
            f'**Unicode Emoji Name:** {unicode_name}\n'
            f'**Code Points:** {code_points}\n'
            f'**Description:** {description}'
        )
        return unicode_info
    else:
        return 'Unicode emoji not found. Make sure it is a valid Unicode emoji character.'

def get_proximity(default, input) -> float:
    default_fp = AllChem.GetMorganFingerprintAsBitVect(default, 2)
    input_fp = AllChem.GetMorganFingerprintAsBitVect(input, 2)
    similarity = DataStructs.FingerprintSimilarity(default_fp, input_fp)
    return similarity

def get_mol(arg):
    try:
        compounds = pcp.get_compounds(arg, 'name')
        compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
        mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
        if mol is None:
            raise ValueError('Invalid SMILES string')
        return mol
    except:
        mol = Chem.MolFromSmiles(arg)
        if mol is None:
            raise ValueError('Invalid SMILES string')
        return mol

def get_molecule_name(molecule) -> str:
    smiles = Chem.MolToSmiles(molecule)
    if smiles:
        compounds = pcp.get_compounds(smiles, 'smiles') #, record_type='3d')
        compound_data = compounds[0].to_dict(properties=['synonyms'])
        if not compounds:
            raise ValueError('No compound found for the given SMILES string')
        return compound_data['synonyms'][0]

def get_sds(query: str):
    query_lower = query.lower()
    try:
        compound = pcp.get_compounds(query, 'name')[0]
        compound_info = compound.to_dict(properties=['iupac_name', 'molecular_formula', 'molecular_weight', 'synonyms'])
        description = compound_info.get('synonyms', ['No description available'])[0]
        pubchem_url = f'https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid}'
        pubmed_description_url = f'https://www.ncbi.nlm.nih.gov/sites/entrez?LinkName=pccompound_pubmed&db=pccompound&cmd=Link&from_uid={compound.cid}'
        description = compound.record['description'][0] if 'description' in compound.record else description
    except (KeyError, IndexError):
        description = 'No detailed description available.'
    embed = discord.Embed(title=description, color=0x00ff99)
    embed.add_field(name='Molecule Name', value=compound_info.get('iupac_name', 'N/A'), inline=False)
    embed.add_field(name='Your Input', value=query, inline=True)
    embed.add_field(name='PubChem URL', value=f'[View PubChem]({pubchem_url})' if pubchem_url else 'No PubChem URL found.', inline=True)
    embed.add_field(name='PubMed URL', value=f'[View PubMed]({pubmed_description_url})' if pubmed_description_url else 'No PubMed URL found.', inline=True)
    return embed

def gsrs(arg):
    chrome_options = Options()
    chrome_options.add_argument('--headless')  # Run headless Chrome (no UI)
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument('--disable-dev-shm-usage')
    driver = webdriver.Chrome(service=Service(executable_path='/home/spawd/.local/bin/chromedriver'), options=chrome_options)
    try:
        search_url = f'https://gsrs.ncats.nih.gov/ginas/app/beta/browse-substance?search={arg}'
        driver.get(search_url)
        driver.implicitly_wait(10)  # Adjust the wait time as needed
        img_element = driver.find_element(By.CSS_SELECTOR, 'body > app-root > app-base > app-substances-browse > div > div.substance-cards > app-substance-summary-card > mat-card > mat-card-title > a')
        if img_element:
            img_src = img_element.get_attribute('href')
            if img_src:
                stripped = img_src.split('/', -1)[-1:]
                link = f'https://gsrs.ncats.nih.gov/api/v1/substances/render({stripped[0]})?format=png&size=512&stereo=true'
                response = requests.get(link)
                image_bytes = response.content
                image = Image.open(BytesIO(image_bytes))
                image = image.convert('RGBA')
                draw = ImageDraw.Draw(image)
                width, height = image.size
                diagonal = math.sqrt(width**2 + height**2)
                font_size = int(diagonal / 15)
                try:
                    font = ImageFont.truetype('Roboto-Regular.ttf', font_size)
                except IOError:
                    font = ImageFont.load_default()
                bbox = draw.textbbox((0, 0), arg, font=font)
                text_width = bbox[2] - bbox[0]
                text_height = bbox[3] - bbox[1]
                text_x = (width - text_width) / 2
                text_y = (height - text_height) / 2
                watermark_image = Image.new('RGBA', image.size, (0, 0, 0, 0))
                watermark_draw = ImageDraw.Draw(watermark_image)
                watermark_draw.text((text_x, text_y), arg, font=font, fill=(255, 255, 255, 64))
                mask = watermark_image.split()[3]
                image.paste(watermark_image, (0, 0), mask)
                return image
            else:
                return 'No src attribute found in the <img> element'
        else:
            return 'No <img> element found with the specified CSS path'
    finally:
        driver.quit()

async def handle_command_error(ctx: commands.Context, error: commands.CommandError):
    log_id = str(uuid.uuid4())
    logger.error(f"Log ID {log_id} - Error in command {ctx.command}: {error}")
    logger.error(f"Full stack trace:\n{traceback.format_exc()}")
    await ctx.send(f"An unexpected error occurred (Log ID: `{log_id}`). Please report this to the support team.")

def increment_version(config: Dict[str, Any]):
    current_version = config['version']
    major, minor, patch = map(int, current_version.split('.'))
    patch += 1
    if patch >= 10:
        patch = 0
        minor += 1
    if minor >= 10:
        minor = 0
        major += 1
    new_version = f'{major}.{minor}.{patch}'
    config['version'] = new_version
    with open(path_config_yaml, 'w') as file:
        yaml.dump(config, file)

def setup_logging(config: Dict[str, Any]) -> None:
    logging_level = config['logging_level'].upper()
    logging.basicConfig(level=getattr(logging, logging_level))
    if not exists(dirname(path_log)):
        makedirs(dirname(path_log))
    file_handler = logging.FileHandler(path_log)
    file_handler.setLevel(logging_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger = logging.getLogger()
    logger.setLevel(logging_level)
    logger.addHandler(file_handler)

def stable_cascade(prompt):
    try:
        client = Client('multimodalart/stable-cascade')
        result = client.predict(
            prompt=prompt,
            negative_prompt='',
            seed=randint(0, 2147483647),
            width=1024,
            height=1024,
            prior_num_inference_steps=20,
            prior_guidance_scale=4,
            decoder_num_inference_steps=10,
            decoder_guidance_scale=0,
            num_images_per_prompt=1,
            api_name="/run",
        )
        return discord.File(result, 'image.webp')
    except ConnectionError as conn_err:
        print(f"Connection error: {conn_err}")
        return "Failed to connect to the server. Please try again later."
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return f"An error occurred: {e}"

def unique_pairs(strings_list):
    pairs = list(itertools.combinations(strings_list, 2))
    sorted_pairs = [sorted(list(pair)) for pair in pairs]
    sorted_pairs_overall = sorted(sorted_pairs)
    return sorted_pairs_overall

def view_difference(molecule1, molecule2):
    mcs_result = rdFMCS.FindMCS([molecule1, molecule2])
    mcs_molecule = Chem.MolFromSmarts(mcs_result.smartsString)
    non_matching_atoms1 = [
        atom.GetIdx() 
        for atom in molecule1.GetAtoms() 
        if atom.GetIdx() not in molecule1.GetSubstructMatch(mcs_molecule)
    ]
    non_matching_atoms2 = [
        atom.GetIdx()
        for atom in molecule2.GetAtoms() 
        if atom.GetIdx() not in molecule2.GetSubstructMatch(mcs_molecule)
    ]
    img = Draw.MolsToGridImage(
        [molecule1, molecule2],
        highlightAtomLists=[non_matching_atoms1, non_matching_atoms2],
        legends=["Molecule 1", "Molecule 2"]
    )
    img_bytes = BytesIO()
    img.save(img_bytes, format='PNG')
    img_bytes.seek(0)  # Go to the beginning of the BytesIO buffer
    
    # Send the image in the Discord channel
    return discord.File(img_bytes, filename='molecules_diff.png')
