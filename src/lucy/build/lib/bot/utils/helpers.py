""" helpers.py
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
"""

from bot.main import Lucy
from datetime import datetime
from os.path import isfile
from discord.ext import commands
from functools import wraps
import asyncio
import discord
import emoji as emoji_lib
import json
import os
import pubchempy as pcp
import random
import requests
import unicodedata

NCBI_REQUEST_DELAY = 0.34  # Delay to ensure no more than 3 requests per second
tasks = {}  # Dictionary to store tasks keyed by (command name, user ID)

def create_unique_pairs_with_indexes(input_bytes_list):
    all_unique_pairs = []
    all_index_pairs = []
    for input_list in input_bytes_list:
        unique_pairs = set()  # Use a set to avoid duplicate pairs
        index_pairs = []  # This will store the index pairs as a list of lists
        for i in range(len(input_list)):
            for j in range(i + 1, len(input_list)):
                pair1 = (input_list[i], input_list[j])
                pair2 = (input_list[j], input_list[i])
                if {0, 1}.issubset({input_list[i], input_list[j]}):
                    continue
                if random.choice([True, False]):
                    if pair1 not in unique_pairs:
                        unique_pairs.add(pair1)
                        index_pairs.append([i, j])  # Append as a list of indexes
                else:
                    if pair2 not in unique_pairs:
                        unique_pairs.add(pair2)
                        index_pairs.append([j, i])  # Append as a list of indexes
        result_list = [list(pair) for pair in unique_pairs]
        all_unique_pairs.append(result_list)
        all_index_pairs.append(index_pairs)
    return all_unique_pairs, all_index_pairs

async def delayed_command(ctx, delay_seconds):
    task_key = (ctx.command.name, ctx.author.id)
    if task_key in tasks and tasks[task_key]:
        tasks[task_key].cancel()
    async def delayed_action():
        await asyncio.sleep(delay_seconds)
        def is_bot_message(message):
            return message.author == self.bot.user
        deleted = await self.purge_messages(ctx, limit, is_bot_message)
    tasks[task_key] = asyncio.create_task(delayed_action())

def fetch_ncbi_organism(tax_id):
    """
    Fetch detailed information about an organism from NCBI Taxonomy database.
    """
#    https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    summary_params = {
        "db": "taxonomy",
        "id": tax_id,
        "retmode": "json"
    }
    response = requests.get(summary_url, params=summary_params)
    data = response.json()
    if response.status_code == 200 and "result" in data and tax_id in data["result"]:
        organism_info = data["result"][tax_id]
        return {
            "tax_id": tax_id,
            "scientific_name": organism_info.get("scientificname", "N/A"),
            "rank": organism_info.get("rank", "N/A"),
            "division": organism_info.get("division", "N/A"),
            "genetic_code": organism_info.get("genetic_code", "N/A"),
            "lineage": organism_info.get("lineage", "N/A"),
            "description": organism_info.get("description", "N/A")
        }
    else:
        return {
            "tax_id": tax_id,
            "scientific_name": "N/A",
            "rank": "N/A",
            "division": "N/A",
            "genetic_code": "N/A",
            "lineage": "N/A",
            "description": "N/A"
        }

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

def get_sds(query: str):
    compound = pcp.get_compounds(query, 'name')[0]
    compound_info = compound.to_dict(properties=['iupac_name', 'molecular_formula', 'molecular_weight', 'synonyms'])
    description = compound_info.get('synonyms', ['No description available'])[0]
    pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{compound.cid}"
#    fda_guide_url = search(query, 'fda')
    msds_message = (
      f"**Molecule Name:** {compound_info.get('iupac_name', 'N/A')}\n"
      f"**Molecular Formula:** {compound_info.get('molecular_formula', 'N/A')}\n"
      f"**Molecular Weight:** {compound_info.get('molecular_weight', 'N/A')} g/mol\n"
      f"**Description:** {description}\n"
      f"**PubChem URL:** [View PubChem]{pubchem_url if pubchem_url else 'No PubChem url found.'}"
 #     f"**FDA:** {fda_guide_url[0]['link'] if fda_guide_url[0]['link'] else 'No prescriber guide found.'}"
    )
    return msds_message

def is_off_peak_hours():
    """
    Determine if the current time is during off-peak hours.
    Off-peak hours: Weekdays from 9:00 PM to 5:00 AM ET and weekends.
    """
    now = datetime.now()
    if now.weekday() >= 5:  # Saturday (5) or Sunday (6)
        return True
    if now.hour > 5 or now.hour <= 21:  # Before 5 AM or after 9 PM
        return True
    return False

def load_config():
    home_dir = os.path.expanduser('~')
    config_path = os.path.join(home_dir, '.config', 'lucy', 'config.json')
    if not os.path.exists(config_path):
        raise FileNotFoundError("Configuration file not found.")
    with open(config_path, 'r') as f:
        return json.load(f)

def off_peak_hours_required():
    """
    Decorator to restrict command execution to off-peak hours.
    """
    def decorator(func):
        @wraps(func)
        async def wrapper(*args, **kwargs):
            self = args[0]  # The first argument is typically `self` or `ctx`
            ctx = args[1]  # The first argument is typically `self` or `ctx`
            if is_off_peak_hours():
                return await func(*args, **kwargs)
            else:
                await ctx.send("This operation is restricted to off-peak hours (9:00 PM - 5:00 AM ET on weekdays, or anytime on weekends). Please try again later.")
        return wrapper
    return decorator

def search_taxonomy(term, max_results=10):
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_params = {
        "db": "taxonomy",
        "term": term,
        "retmode": "json",
        "retmax": max_results
    }
    search_response = requests.get(search_url, params=search_params)
    search_data = search_response.json()
    if not search_data["esearchresult"]["idlist"]:
        return []
    tax_ids = search_data["esearchresult"]["idlist"]
    return [{"tax_id": tax_id} for tax_id in tax_ids]
