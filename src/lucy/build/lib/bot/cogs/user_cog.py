''' user_cog.py
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
from discord.ext import commands
from PIL import Image, ImageFont, ImageDraw
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, DataStructs, Draw
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.DataStructs import FingerprintSimilarity, TanimotoSimilarity
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
from rdkit.Chem.Draw import IPythonConsole

import colorsys
import discord
from io import BytesIO
import io
import math
import os
import PIL
import pubchempy as pcp
import rdkit
import traceback
import wikipediaapi

from bot.utils.helpers import delayed_command
from bot.utils.helpers import fetch_ncbi_organism
from bot.utils.helpers import get_emoji
from bot.utils.helpers import get_sds
from bot.utils.helpers import load_config
from bot.utils.helpers import off_peak_hours_required
from bot.utils.helpers import search_taxonomy
from discord.ext import commands
from functools import wraps
from googletrans import Translator, LANGUAGES
from manim import *
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager
from PIL import Image, ImageFont, ImageDraw

from googleapiclient.discovery import build
import random
import discord
import json
import io
import math
import os
import PIL
import pubchempy as pcp
import requests
import tempfile
import traceback
from bot.utils.helpers import create_unique_pairs_with_indexes


class UserCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.translator = Translator()
        self.user_translation_preferences = {}
#        self.user_id = '154749533429956608'
        self.api_url = 'https://api.psychonautwiki.org/'
        self.tags = {}
        self.guild = self.bot.get_guild(self.bot.testing_guild_id)
        self.channels = bot.get_all_channels()
        self.headers = {
            "User-Agent": "DiscordBot/1.0 (https://github.com/brandongrahamcobb/Documents; brandongrahamcobb@icloud.com)"
        }

    @commands.Cog.listener()
    async def on_command_error(self, ctx: commands.Context, error):
        if isinstance(error, commands.CommandOnCooldown):
            await ctx.send(f'Command is on cooldown. Try again in {error.retry_after:.2f} seconds.')
        else:
            raise error

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message):
        if message.author.bot:
            return
        prefs = self.user_translation_preferences.get(message.author.id)
        if prefs:
            target_lang_code, source_lang_code = prefs
            try:
                translated = self.translator.translate(message.content, dest=target_lang_code, src=source_lang_code)
                await message.channel.send(f'{translated.text}')
            except Exception as e:
                await message.channel.send(f'Error translating message: {e}')

    def add_watermark(self, image: BytesIO, watermark_text: str = 'Unknown') -> BytesIO:
        RGB_image = Image.open(image)
        RGBA_image = RGB_image.convert('RGBA')
        draw = ImageDraw.Draw(RGBA_image)
        width, height = RGBA_image.size
        diagonal = math.sqrt(width**2 + height**2)
        font_size = int(diagonal / 15)
        try:
            font = ImageFont.truetype("/home/spawd/Documents/src/CobbBrandonGraham_lucy_220824/resources/Roboto-Regular.ttf", font_size)  # Replace with the path to your font file
        except IOError:
            font = ImageFont.load_default()
        while True:
            bbox = draw.textbbox((0, 0), watermark_text, font=font)
            text_width = bbox[2] - bbox[0]
            if text_width <= 512 or font_size <= 1:
                break
            font_size -= 1
            font = ImageFont.truetype("/home/spawd/Documents/src/CobbBrandonGraham_lucy_220824/resources/Roboto-Regular.ttf", font_size)  # Replace with the path to your font file
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

    def adjust_hue_and_saturation(self, image, hue_shift, saturation_shift) -> BytesIO:
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

    def combine(self, bytes1: BytesIO, bytes2: BytesIO) -> BytesIO:
        img1 = Image.open(bytes1) #inv_img1 = self.invert_colors(bytes1)
        #inv_wimg1 = Image.open(self.add_watermark(inv_img1, watermark_text=resolved_name))                          ###############
        #inv_img2 = self.invert_colors(bytes2)
        img2 = Image.open(bytes2) #inv_img1 = self.invert_colors(bytes1)
       # inv_wimg2 = Image.open(self.add_watermark(inv_img2, watermark_text=resolved_name))
        widths, heights = zip(*(img.size for img in [img1, img2]))
        total_width = sum(widths)
        max_height = max(heights)
        combined_img = Image.new('RGB', (total_width, max_height))
        x_offset = 0
        for img in [img1, img2]:
            combined_img.paste(img, (x_offset, 0))
            x_offset += img.width
        inverted_image = self.invert_colors(combined_img)
        output = self.adjust_hue_and_saturation(inverted_image, hue_shift=-180, saturation_shift=160)
        output.seek(0)
        return output

    def draw_fingerprint(self, queries) -> BytesIO:
        molecules = []
        for smiles in queries:
            mol = Chem.MolFromSmiles(self.get_smiles(smiles))
            molecules.append(mol)
        compounds = pcp.get_compounds(self.get_smiles(queries[0]), 'smiles')
        if compounds[0].synonyms == None:
            resolved_name = '~LUCY~'
        else:
            resolved_name = compounds[0].synonyms[0]
        d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
        d2d.prepareMolsBeforeDrawing = False
        Options = d2d.drawOptions()
        Options.prepareMolsBeforeDrawing = False
        Options.includeMetadata = False
        d2d.SetDrawOptions(Options)
        mol1 = rdMolDraw2D.PrepareMolForDrawing(molecules[0], kekulize=True)
        mol1.UpdatePropertyCache(False)
        mol2 = rdMolDraw2D.PrepareMolForDrawing(molecules[1], kekulize=True)
        mol2.UpdatePropertyCache(False)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(molecules[0], 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(molecules[1], 2, nBits=2048)
        similarity = TanimotoSimilarity(fp1, fp2)
        fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol1, mol2, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d, drawingOptions=Options, legend=resolved_name)
        d2d.FinishDrawing()
        drawing = d2d.GetDrawingText()
        output = BytesIO(drawing) #output = self.add_watermark(BytesIO(drawing), watermark_text=resolved_name)
        return output

    def draw_watermarked_molecule(self, arg: str) -> BytesIO:
        smiles = self.get_smiles(arg)
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                compounds = pcp.get_compounds(smiles, 'smiles')
                if compounds[0].synonyms == None:
                    resolved_name = '~LUCY~'
                else:
                    resolved_name = compounds[0].synonyms[0]
                d2d = rdMolDraw2D.MolDraw2DCairo(512, 512)
                Options = d2d.drawOptions()
                d2d.SetDrawOptions(Options)
                rdMolDraw2D.SetDarkMode(Options)
                mol1 = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=True)
                mol1.UpdatePropertyCache(False)
                d2d.DrawMolecule(mol1)
                d2d.FinishDrawing()
                drawing = d2d.GetDrawingText()
                output = self.add_watermark(BytesIO(drawing), watermark_text=resolved_name)
                return output
            else:
                return None
        else:
            return None

    def fetch_taxonomy(self, taxon_name):
        sparql_arg = f"""
        SELECT ?taxon ?taxonLabel ?image ?description WHERE {{
          wd:{taxon_name} wdt:P171 ?parent .
          ?taxon wdt:P171 wd:{taxon_name} .
          OPTIONAL {{ ?taxon wdt:P18 ?image. }}
          OPTIONAL {{ ?taxon schema:description ?description. FILTER (lang(?description) = "en") }}
          SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
        }}
        """
        response = requests.get('https://arg.wikidata.org/sparql', headers=self.headers, params={"arg": sparql_arg, "format": "json"})
        return response.json() if response.status_code == 200 else None

    def get_language_code(self, language_name):
        language_name = language_name.lower()
        for lang_code, lang_name in LANGUAGES.items():
            if lang_name.lower() == language_name:
                return lang_code
        return None

    def get_smiles(self, arg: str) -> str:
        mol = Chem.MolFromSmiles(arg) #or Chem.MolFromInchl(arg)
        if mol:
            return arg
        else:
            results = pcp.get_compounds(arg, 'name')
            if results:
                return results[0].isomeric_smiles
            else:
                return None

    def gsrs(self, arg):
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
                    image = Image.open(io.BytesIO(image_bytes))
                    image = image.convert('RGBA')
                    draw = ImageDraw.Draw(image)
                    width, height = image.size
                    diagonal = math.sqrt(width**2 + height**2)
                    font_size = int(diagonal / 15)  # Adjust this factor as needed
                    try:
                        font = ImageFont.truetype('arial.ttf', font_size)
                    except IOError:
                        font = ImageFont.load_default(font_size)
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

    def invert_colors(self, image):
        return Image.eval(image, lambda x: 255 - x)

    async def on_cooldown(self, ctx, command_name: str):
        command = self.bot.get_command(command_name)
        if not command:
            await ctx.send(f"No command found with the name `{command_name}`.")
            return True
        bucket = command._buckets.get_bucket(ctx)
        retry_after = bucket.update_rate_limit()
        if retry_after:
            await ctx.send(f"Command `{command_name}` is on cooldown. Try again in {retry_after:.2f} seconds.")
            return True
        return False

    def proximity(self, default: str, input: str, flag: bool = True):
        similarity_list = []
        input_compounds = pcp.get_compounds(self.get_smiles(input), 'smiles')
        input_smiles = input_compounds[0].isomeric_smiles
        input_mol = Chem.MolFromSmiles(input_smiles)
        input_fp = AllChem.GetMorganFingerprintAsBitVect(input_mol, 2)

        default_compounds = pcp.get_compounds(self.get_smiles(default), 'smiles')
        default_smiles = default_compounds[0].isomeric_smiles
        default_mol = Chem.MolFromSmiles(default_smiles)
        default_fp = AllChem.GetMorganFingerprintAsBitVect(default_mol, 2)
        default_tanimoto = DataStructs.FingerprintSimilarity(default_fp, input_fp)

        if flag:
            LSD = 'CCN(CC)C(=O)[C@H]1CN([C@@H]2CC3=CNC4=CC=CC(=C34)C2=C1)C'
            lsd_mol = Chem.MolFromSmiles(LSD)
            lsd_fp = AllChem.GetMorganFingerprintAsBitVect(lsd_mol, 2)
            lsd_tanimoto = DataStructs.FingerprintSimilarity(lsd_fp, input_fp)
            is_closer_to_default = default_tanimoto > lsd_tanimoto

        similarity = TanimotoSimilarity(default_fp, input_fp)
        similarity_list.append(similarity)
#            unique_similarity_list = create_unique_psimilarity_list)
        return similarity_list

#    @commands.command(name='animate', hidden=True)
#    async def animate(self, ctx, *, expression: str):
#        await ctx.send(f"Generating animation for: {expression}")
#        class MyScene(Scene):
#            def construct(self):
#                formula = MathTex(expression)
#                self.play(Write(formula))
#                self.wait(2)
#        with tempfile.TemporaryDirectory() as tempdir:
#            config.output_file = os.path.join(tempdir, "animation.mp4")
#            scene = MyScene()
#            scene.render()
#            video_path = config.output_file
#            with open(video_path, 'rb') as video_file:
#                video_bytes = video_file.read()
#            video_io = io.BytesIO(video_bytes)
#            video_io.seek(0)
#            await ctx.send(file=discord.File(video_io, filename="animation.mp4"))


    @commands.command(name="covid", help="Get the latest COVID-19 statistics globally or for a specific country")
    async def covid(self, ctx, country: str = 'global'):
        """Fetches COVID-19 statistics from the disease.sh API."""
        if country.lower() == 'global':
            url = 'https://disease.sh/v3/covid-19/all'
        else:
            url = f'https://disease.sh/v3/covid-19/countries/{country}'

        response = requests.get(url)
        data = response.json()

        if response.status_code == 200:
            embed = discord.Embed(title=f"COVID-19 Stats for {country.title()}", color=discord.Color.blue())
            embed.add_field(name="Cases", value=f"{data['cases']:,}", inline=True)
            embed.add_field(name="Deaths", value=f"{data['deaths']:,}", inline=True)
            embed.add_field(name="Recovered", value=f"{data['recovered']:,}", inline=True)
            embed.add_field(name="Active", value=f"{data['active']:,}", inline=True)
            embed.add_field(name="Critical", value=f"{data['critical']:,}", inline=True)
            embed.add_field(name="Today's Cases", value=f"{data['todayCases']:,}", inline=True)
            embed.add_field(name="Today's Deaths", value=f"{data['todayDeaths']:,}", inline=True)
            await ctx.send(embed=embed)
        else:
            await ctx.send("Failed to retrieve data. Please try again later or check the country name.")

    @commands.cooldown(rate=3, per=1.0, type=commands.BucketType.channel)
    @commands.command(name='draw')
    @off_peak_hours_required()
    async def draw(self, ctx, *args) -> None:
        if await self.on_cooldown(ctx, 'draw'):
            return
        async with ctx.typing():
            molecules = []
            for arg in args:
                image = self.draw_watermarked_molecule(arg)
                if image:
                    molecules.append(image)
            if len(molecules) == 1:
                for i, image in enumerate(molecules):
                  #  embed = discord.Embed()
         #           for channel in self.channels:
          #              if channel.name == 'vegan':
                 #   await ctx.send(file=discord.File(image, f'molecule_{i+1}.png'))
                    link = await ctx.send(file=discord.File(image, f'molecule_{i+1}.png'))
            #        await ctx.send(f'I produced {args[i]} over on my test guild!' + link.jump_url)
                            #embed.add_field(name=f'I produced {args[i]} over on my test guild!', value=link.jump_url)
                   #  await ctx.send('Here is a link to my server: https://discord.gg/dNfUn8MeYN')
                   # await ctx.send(embed=embed)
            elif len(molecules) >= 2:
                all_unique_pairs, all_index_pairs = create_unique_pairs_with_indexes(args)
                pair_of_names = []
                for pair in range(len(all_unique_pairs)):
                    pair_of_names.append(args[all_unique_pairs[pair][0])
                    pair_of_names.append(args[all_unique_pairs[pair][1])
                    proximity = self.proximity(default=pair_of_names[0], input=pair_of_names[1], flag=True)
                    await ctx.send(f'The Tanimoto Similarity of {pair_of_names[0]} and {pair_of_names[1]} is {proximity:.2f}.')
                    fingerprints = []
                    fingerprints.append(self.draw_fingerprint(list(reversed(pair_of_names))))
                    fingerprints.append(self.draw_fingerprint(pair_of_names))
                    combined_image = self.combine(fingerprints[0], fingerprints[1])
                    await ctx.send(file=discord.File(combined_image, f'molecule_comparison.png'))
                    pair_of_names = []
            else:
                await ctx.send('No valid molecules found.')

#    @commands.command(name='get', hidden=True)
#    @commands.is_owner()
#    async def get(self, ctx, *, arg: str):
#        results = search(arg, 'msds')
#        embed = discord.Embed()
#        for result in results:
#            embed.add_field(name=result['title'], value=result['link'], inline=False)
#        await ctx.send(embed=embed)

    @commands.command(name='get', hidden=True)
    @commands.is_owner()
    async def random_plant(self, ctx: commands.Context, search_arg: str = 'cannabis cola'):
        google_json_key = (load_config())['api_keys']['api_key_1']
        cse_id = 'e1753e22095c74f10'  # Using your 'web' CSE ID
        all_images = []
        start_index = 1
        max_requests = 5  # Adjust this to control how many requests you want to make
        service = build('customsearch', 'v1', developerKey=google_json_key)
        for _ in range(max_requests):
            # Perform a Google Image search with pagination
            res = service.cse().list(q=search_arg, cx=cse_id, searchType='image', num=10, start=start_index).execute()
            images = res.get('items', [])
            if not images:
                break
            all_images.extend(images)
            start_index += 10  # Move to the next page of results
        valid_images = []
        for image in all_images:
            image_url = image['link']
            try:
                response = requests.get(image_url)
                img = Image.open(BytesIO(response.content))
                width, height = img.size
                if width >= 1920 and height >= 1080 and width / height == 16 / 9:
                    valid_images.append(image_url)
            except Exception as e:
                print(f"Error processing image {image_url}: {e}")
        if not valid_images:
            await ctx.send("No suitable images found.")
            return
        random_image = random.choice(valid_images)
        await ctx.send(random_image)

    @commands.command(name='info', hidden=True)
    async def info(self, ctx: commands.Context, argument):
        if argument is None:
            await ctx.send('The command is used !info <MOLECULE> or !info :emoji:')
            return
        if isinstance(argument, discord.Emoji):
            try:
                await ctx.send(get_emoji(argument))
            except Exception as e:
                await ctx.send(e)
        else:
            try:
                await ctx.send(get_sds(argument))
            except Exception as e:
                await ctx.send(e)

    @commands.command(name='languages', hidden=True)
    @commands.is_owner()
    async def languages(self, ctx):
        supported_languages = ', '.join(LANGUAGES.values())
        await ctx.send(f'Supported languages are:\n{supported_languages}')

    @commands.command(name='logp', hidden=True)
    @commands.is_owner()
    async def logp(self, ctx, *args):
        try:
            for arg in args:
                compounds = pcp.get_compounds(arg, 'name')
                compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
                mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
                log_p = Crippen.MolLogP(mol)
                await ctx.send(f'Your octanol:water coefficient is: {log_p}')
        except Exception as e:
            await ctx.send(f'Error fetching data: {e}')

    @commands.command(name='color')
    async def role_or_pdf(self, ctx, role_or_file: discord.Role = None):
        if isinstance(role_or_file, discord.Role):
            # Get the role's color
            color = role_or_file.color
            await ctx.send(f"The RGB color of the role `{role_or_file.name}` is: {color.to_rgb()}")
        elif ctx.message.attachments:
            attachment = ctx.message.attachments[0]
            if attachment.filename.lower().endswith('.pdf'):
                # Download the PDF file
                pdf_data = await attachment.read()

                # Convert first page of PDF to PNG
                with io.BytesIO(pdf_data) as pdf_io:
                    pdf_reader = PdfFileReader(pdf_io)
                    first_page = pdf_reader.getPage(0)

                    # Convert the first page to an image using PIL
                    # For demonstration, I'll use a placeholder image
                    image = Image.new('RGB', (128, 128), color='white')
                    draw = ImageDraw.Draw(image)
                    draw.text((32, 32), "PDF Page", fill='black')

                    # Resize to 128x128
                    image = image.resize((128, 128))

                    # Save the image to a bytes object
                    with io.BytesIO() as image_binary:
                        image.save(image_binary, 'PNG')
                        image_binary.seek(0)

                        # Upload the image as an emoji
                        emoji = await ctx.guild.create_custom_emoji(name='pdf_page', image=image_binary.read())
                        await ctx.send(f"Created emoji: <:{emoji.name}:{emoji.id}>")
            else:
                await ctx.send("Please attach a PDF file.")
        else:
            await ctx.send("Please provide a valid role or attach a PDF file.")

    @commands.command(name="s", hidden=True)
    #@commands.is_owner()
    async def s(self, ctx: commands.Context, *molecules):
    #    if True in booleans:
            #command = self.bot.get_command('info')
            #await ctx.invoke(command, argument=molecule)
   #         await ctx.send(f'What does Spawd think of {molecule}?\n\N{OK HAND SIGN}')
  #      else:
 #           await ctx.send(f'what does Spawd think of {molecule}?\n\N{REVERSED HAND WITH MIDDLE FINGER EXTENDED}')
#        for floats

#        ativan = 'C1=CC=C(C(=C1)C2=NC(C(=O)NC3=C2C=C(C=C3)Cl)O)Cl'
        nicotine = 'CN1CCC[C@H]1C2=CN=CC=C2'
        B12 = '[CH3-].CC1=CC2=C(C=C1C)N(C=N2)[C@@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)([O-])OC(C)CNC(=O)CC[C@@]4([C@H]([C@@H]5[C@]6([C@@]([C@@H](/C(=C(/C7=N/C(=C\C8=N/C(=C(\C4=N5)/C)/[C@H](C8(C)C)CCC(=O)N)/[C@H]([C@]7(C)CC(=O)N)CCC(=O)N)\C)/[N-]6)CCC(=O)N)(C)CC(=O)N)C)CC(=O)N)C)O.[Co+3]'
        caffeine = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
        catapres = 'C1CN=C(N1)NC2=C(C=CC=C2Cl)Cl'
        CBD = 'CCCCCC1=CC(=C(C(=C1)O)[C@@H]2C=C(CC[C@H]2C(=C)C)C)O'
        creatine = 'CN(CC(=O)O)C(=N)N'
        D3 = 'C[C@H](CCCC(C)C)[C@H]1CC[C@@H]\2[C@@]1(CCC/C2=C\C=C/3\C[C@H](CCC3=C)O)C'
        DHA = 'CC/C=C\C/C=C\C/C=C\C/C=C\C/C=C\C/C=C\CCC(=O)O'
        DPA = 'CC/C=C\C/C=C\C/C=C\C/C=C\C/C=C\CCCCCC(=O)O'
        luvox = 'COCCCC/C(=N\OCCN)/C1=CC=C(C=C1)C(F)(F)F.C(=C\C(=O)O)\C(=O)O'
        SAME = 'C[S+](CC[C@@H](C(=O)[O-])N)C[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)O'
        speed = 'C[C@@H](CC1=CC=CC=C1)N'
        THCA = 'CCCCCC1=CC2=C([C@@H]3C=C(CC[C@H]3C(O2)(C)C)C)C(=C1C(=O)O)O'
        THCA = 'CCCCCC1=CC2=C([C@@H]3C=C(CC[C@H]3C(O2)(C)C)C)C(=C1C(=O)O)O'
        triamcinolone_acetonide = 'C[C@]12C[C@@H]([C@]3([C@H]([C@@H]1C[C@@H]4[C@]2(OC(O4)(C)C)C(=O)CO)CCC5=CC(=O)C=C[C@@]53C)F)O'
        vraylar = 'CN(C)C(=O)NC1CCC(CC1)CCN2CCN(CC2)C3=C(C(=CC=C3)Cl)Cl'
        ativan = 'C1=CC=C(C(=C1)C2=NC(C(=O)NC3=C2C=C(C=C3)Cl)O)Cl'
#            list.append(self.proximity(ctx, molecule, B12))
#        sorted_a = create_unique_pairs(a)
 #           list.append(self.proximity(ctx, molecule, caffeine))
 #           list.append(self.proximity(ctx, molecule, catapres))
#            list.append(self.proximity(ctx, molecule, CBD))
   #         list.append(self.proximity(ctx, molecule, creatine))
  #          list.append(self.proximity(ctx, molecule, D3))
     #       list.append(self.proximity(ctx, molecule, DHA))
      #      list.append(self.proximity(ctx, molecule, DPA))
            #a, b = self.proximity(ctx, molecule, luvox))
            #a, b = self.proximity(ctx, molecule, THCA))
            #a, b = list.append(self.proximity(ctx, molecule, vraylar))
            #list.append(a)
#            list.append(self.proximity(ctx, molecule, nicotine))
          #  list.append(self.proximity(ctx, molecule, ativan))
        #    list.append(self.proximity(ctx, molecule, SAME))
#            list.append(self.proximity(ctx, molecule, speed))
      #      list.append(self.proximity(ctx, molecule, triamcinolone_acetonide))

#    @commands.command(name='tag', hidden=True)
#    @commands.is_owner()
#    async def tag(self, ctx, identifier: str):
#        if ctx.message.attachments:
#            attachment = ctx.message.attachments[0]
#            file_path = f"resources/{identifier}_{attachment.filename}"
#            await attachment.save(file_path)
#            self.tags[identifier] = file_path
#            await ctx.send(f"Image tagged as '{identifier}' successfully!")
#        else:
#            file_path = self.tags.get(identifier)
#            if file_path and os.path.exists(file_path):
#                await ctx.send(file=discord.File(file_path))
#            else:
#                await ctx.send(f"No image found with the identifier '{identifier}'.")


    @commands.command(name='search', hidden=True)
    async def search(self, ctx, *, arg=None):
        try:
            try:
                watermarked_image = self.gsrs(arg)
                with io.BytesIO() as image_binary:
                    watermarked_image.save(image_binary, format='PNG')
                    image_binary.seek(0)
                    await ctx.send(file=discord.File(fp=image_binary, filename='watermarked_image.png'))
            except:
                await ctx.send(traceback.format_exc())
        except:
            if arg is None:
                return
            await ctx.send(f'{arg} is an unknown molecule.')

    @commands.command(name='smiles', hidden=True)
    async def smiles(self, ctx: commands.Context, *, chemical_name: str):
        compounds = pcp.get_compounds(chemical_name, 'name')
        compound = compounds[0]
        isomeric_smiles = compound.isomeric_smiles
        await ctx.send(f'The isomeric SMILES for {chemical_name} is: {isomeric_smiles}')

    @commands.command(name="substance", hidden=True)
    @commands.is_owner()
    async def substance(self, ctx, *, substance_name: str):
        """
        Fetch and display concise information about a substance from PsychonautWiki.
        """
        arg = {
            "arg": f"""
            {{
                substances(arg: "{substance_name}") {{
                    name
                    roas {{
                        name
                        dose {{
                            units
                            common {{ min max }}
                        }}
                        duration {{
                            peak {{ min max units }}
                        }}
                    }}
                    effects {{
                        name
                        url
                    }}
                }}
            }}
            """
        }
        try:
            response = requests.post(self.api_url, json=arg)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            await ctx.send(f"Error fetching data: {e}")
            return
        data = response.json()
        substances = data.get('data', {}).get('substances', [])
        if not substances:
            await ctx.send(f"No data found for substance: {substance_name}")
            return
        substance_info = substances[0]
        embed = discord.Embed(
            title=substance_info['name'],
            color=discord.Color.blue()
        )
        for roa in substance_info['roas']:
            roa_info = (
                f"**Common Dose:** {roa['dose']['common']['min']}-{roa['dose']['common']['max']} {roa['dose']['units']}\n"
                f"**Peak Duration:** {roa['duration']['peak']['min']}-{roa['duration']['peak']['max']} {roa['duration']['peak']['units']}"
            )
            embed.add_field(name=f"ROA: {roa['name']}", value=roa_info, inline=False)
        max_effects_to_display = 5  # Limit the number of effects to display
        effects_info = ""
        effects_count = 0
        for effect in substance_info['effects']:
            effect_text = f"[{effect['name']}]({effect['url']})\n"
            if len(effects_info) + len(effect_text) > 1024 or effects_count >= max_effects_to_display:
                embed.add_field(name="Effects", value=effects_info, inline=False)
                effects_info = effect_text  # Start a new field
                effects_count += 1
            else:
                effects_info += effect_text
            if effects_count >= max_effects_to_display:
                break
        if effects_info:
            embed.add_field(name="Effects", value=effects_info, inline=False)
        if len(embed) > 6000:
            await ctx.send("The information is too long to display in a single message.")
            return
        await ctx.send(embed=embed)

    @commands.command(name='translate', hidden=True)
#    @commands.is_owner()
    async def translate(self, ctx, toggle: str, target_lang: str = 'english', source_lang: str = 'auto'):
        if toggle.lower() == 'on':
            target_lang_code = self.get_language_code(target_lang)
            source_lang_code = self.get_language_code(source_lang)
            if target_lang_code is None or source_lang_code is None:
                await ctx.send(f'{ctx.author.mention}, please specify valid language names.')
                return
            self.user_translation_preferences[ctx.author.id] = (target_lang_code, source_lang_code)
            await ctx.send(f'{ctx.author.mention}, translation enabled from {source_lang} to {target_lang}.')
        elif toggle.lower() == 'off':
            self.user_translation_preferences[ctx.author.id] = None
            await ctx.send(f'{ctx.author.mention}, translation disabled.')
        else:
            await ctx.send(f'{ctx.author.mention}, please specify `on` or `off`.')

    def fetch_taxonomy(self, taxon_name):
        sparql_arg = f"""
        SELECT ?taxon ?taxonLabel ?image ?description WHERE {{
          wd:{taxon_name} wdt:P171 ?parent .
          ?taxon wdt:P171 wd:{taxon_name} .
          OPTIONAL {{ ?taxon wdt:P18 ?image. }}
          OPTIONAL {{ ?taxon schema:description ?description. FILTER (lang(?description) = "en") }}
          SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
        }}
        """
        return arg_wikidata(sparql_arg)

    @commands.command(name='bio')
    async def taxonomy(self, ctx, *, taxon_name: str):
        data = self.fetch_taxonomy(taxon_name)
        if not data:
            await ctx.send("Could not retrieve data. Please check the taxon name.")
            return
        bindings = data.get('results', {}).get('bindings', [])
        if not bindings:
            await ctx.send("No data found for the given taxon.")
            return
        children = []
        for entry in bindings:
            taxon_label = entry.get('taxonLabel', {}).get('value', 'Unknown')
            taxon_id = entry.get('taxon', {}).get('value', '').split('/')[-1]
            image_url = entry.get('image', {}).get('value', None)
            description = entry.get('description', {}).get('value', 'No description available.')
            embed = discord.Embed(title=taxon_label, description=description)
            if image_url:
                embed.set_image(url=image_url)
            await ctx.send(embed=embed)
            children.append((taxon_label, taxon_id))
        if children:
            child_options = "\n".join([f"{index + 1}. {child[0]}" for index, child in enumerate(children)])
            await ctx.send(f"Choose a taxon to explore further:\n{child_options}")
            def check(m):
                return m.author == ctx.author and m.channel == ctx.channel and m.content.isdigit() and 1 <= int(m.content) <= len(children)
            try:
                msg = await self.bot.wait_for('message', check=check, timeout=60.0)
                selected_index = int(msg.content) - 1
                selected_taxon_id = children[selected_index][1]
                await self.taxonomy(ctx, taxon_name=selected_taxon_id)
            except (ValueError, TimeoutError):
                await ctx.send("Invalid selection or timeout. Please try the command again.")
        else:
            await ctx.send(f"No further children available for **{taxon_name}**.")

async def setup(bot: commands.Bot):
    await bot.add_cog(UserCog(bot))
