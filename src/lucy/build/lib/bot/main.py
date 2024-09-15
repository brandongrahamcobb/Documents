
#    Lucy, a discord.py bot, is an open-source package containing four cogs and a helper file.
#    Copyright (C) 2024  Cobb, Brandon Graham

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public license as published by
#    the Free Software Foundation, either version 3 of the license, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public license for more details.

#    You should have received a copy of the GNU General Public license
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

from discord.ext import commands
from typing import List, Optional, Dict, Any
from os import getenv, makedirs
from os.path import abspath, dirname, expanduser, exists, isfile, join
#from bot.cogs.tag_cog import TagCog

import aiomysql
import asyncio
import datetime as dt
import discord
import json
import logging
import logging.handlers
import os

global logger

def prompt_for_values(prompt: str, default_value: str) -> str:
    value = input(f'{prompt} [{default_value}]: ')
    return value if value else default_value

class Lucy(commands.Bot):

    _config = None

    def __init__(
        self,
        *args,
        initial_extensions: List[str],
        db_pool: aiomysql.Pool,
        testing_guild_id: Optional[int] = None,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.db_pool = db_pool
        self.testing_guild_id = testing_guild_id
        self.initial_extensions = initial_extensions

    @classmethod
    def _get_config(cls) -> Dict[str, Any]:
        current_date = dt.datetime.now().strftime('%d%m%y')
        home = expanduser('~')
        if os.name == 'nt':
            base_path = join(home, 'Program Files', 'lucy-package', 'src', f'CobbBrandonGraham_lucy_{current_date}')
            license = join(base_path, 'license')  # ON PURPOSE 'license' NOT 'path_dest_license'
            path_dest_json = join(getenv('APPDATA'), '.config', 'CobbBrandonGraham', 'config.json')
            path_dest_log = join(getenv('APPDATA'), '.log', 'discord.log')
            cogs_path = join(base_path, 'bot', 'cogs')
        else:
            base_path = join('usr', 'share', 'lucy')
            license = join(base_path, 'license')  # ON PURPOSE 'license' NOT 'path_dest_license'
            path_dest_json = join(home, '.config', 'CobbBrandonGraham', 'config.json')
            path_dest_log = join(home, '.log', 'discord.log')
            cogs_path = join(base_path, 'bot', 'cogs')
        if cls._config is None:
            if isfile(path_dest_json):
                with open(path_dest_json, 'r') as file:
                    data = json.load(file)
                data['cogs'] = prompt_for_values('Enter the cogs.', data.get('cogs', [
                    'bot.cogs.admin_cog',
                    'bot.cogs.game_cog',
                    'bot.cogs.my_cog',
                    'bot.cogs.tag_cog',
                    'bot.cogs.user_cog'
                ]))
                data['command_prefix'] = prompt_for_values('Enter the command prefix', data.get('command_prefix', '!'))
                data['database_url'] = prompt_for_values('Enter the database URL', data.get('database_url', 'brandongcobb.com'))
                data['intents'] = prompt_for_values('Enter the intents', data.get('intents', 'discord.Intents.all()'))
                data['logging_level'] = prompt_for_values('Enter the logging level', data.get('logging_level', 'INFO'))
                data['owner_id'] = prompt_for_values('Enter the owner ID', data.get('owner_id', '154749533429956608'))
                data['testing_guild_id'] = prompt_for_values('Enter the testing guild ID', data.get('testing_guild_id', '1217326055111655507'))
                data['token'] = prompt_for_values('Enter the bot token', data.get('token', ''))
                data['user_agent'] = prompt_for_values('Enter the User-Agent header', data.get('user_agent', 'Lucy'))
                data['version'] = prompt_for_values('Enter the bot version', data.get('version', '1.0.0'))
                data['xp_rate'] = prompt_for_values('Enter the XP rate', data.get('xp_rate', '1'))
                data['api_keys'] = data.get('api_keys', {})
                for i in range(1, 21):
                    key = f'api_key_{i}'
                    current_key = data['api_keys'].get(key, '')
                    data['api_keys'][key] = prompt_for_values(f'Enter API key {i}', current_key)
            else:
                makedirs(dirname(path_dest_json), exist_ok=True)
                data = {
                    'api_keys': {},
                    'cogs': prompt_for_values('Enter the cogs.', [
                        'bot.cogs.admin_cog',
                        'bot.cogs.game_cog',
                        'bot.cogs.my_cog',
                        'bot.cogs.tag_cog',
                        'bot.cogs.user_cog'
                    ]),
                    'command_prefix': prompt_for_values('Enter the command prefix', '!'),
                    'database_url': prompt_for_values('Enter the database URL', 'brandongcobb.com'),
                    'intents': prompt_for_values('Enter the intents', 'discord.Intents.all()'),
                    'logging_level': prompt_for_values('Enter the logging level', 'INFO'),
                    'owner_id': prompt_for_values('Enter the your user ID', '154749533429956608'),
                    'testing_guild_id': prompt_for_values('Enter the testing guild ID', '1217326055111655507'),
                    'token': prompt_for_values('Enter the bot token', ''),
                    'version': prompt_for_values('Enter the bot version', '1.0.0'),
                    'user_agent': prompt_for_values('Enter the User-Agent header', 'Lucy'),
                    'xp_rate': prompt_for_values('Enter the XP rate', '1'),
                }
                for i in range(1, 21):
                    data['api_keys'][f'api_key_{i}'] = prompt_for_values(f'Enter API key {i}', '')
            # Increment version
            current_version = data['version']
            major, minor, patch = map(int, current_version.split('.'))
            patch += 1
            if patch >= 10:
                patch = 0
                minor += 1
            if minor >= 10:
                minor = 0
                major += 1
            new_version = f'{major}.{minor}.{patch}'
            data['version'] = new_version
            with open(path_dest_json, 'w') as file:
                json.dump(data, file, indent=4)
            cls._config = data
        return cls._config

    def setup_logging(self, logging_level: str, path_dest_log: str) -> None:
        logging_level = logging_level.upper()
        logging.basicConfig(level=getattr(logging, logging_level))
        if not exists(dirname(path_dest_log)):
            makedirs(dirname(path_dest_log))
        file_handler = logging.FileHandler(path_dest_log)
        file_handler.setLevel(logging_level)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger = logging.getLogger()
        logger.setLevel(logging_level)
        logger.addHandler(file_handler)

    async def setup_hook(self) -> None:
        for cog in self.initial_extensions:
            await self.load_extension(cog)
        if self.testing_guild_id:
            guild = discord.Object(self.testing_guild_id)
            self.tree.copy_global_to(guild=guild)
            await self.tree.sync(guild=guild)

    async def load_config(self) -> Dict[str, Any]:
        return self._get_config()

async def main():
    # Initialize the database connection pool
    db_pool = await aiomysql.create_pool(
        user='lucy',
        password='lucy',  # Make sure to set your actual password
        host='localhost',
        port=3306,
        db='lucy_db',
        autocommit=True,
        loop=asyncio.get_event_loop()
    )
    # Create and run the bot
    async with Lucy(
        command_prefix=commands.when_mentioned_or('!'),
        db_pool=db_pool,
        intents=discord.Intents.all(),
        initial_extensions=[
            'bot.cogs.admin_cog',
            'bot.cogs.game_cog',
            'bot.cogs.my_cog',
            'bot.cogs.tag_cog',
            'bot.cogs.user_cog'
        ],
        testing_guild_id=int('1217326055111655507')  # Convert to int
    ) as bot:
        config = await bot.load_config()
        current_date = dt.datetime.now().strftime('%d%m%y')
        home = expanduser('~')
        if os.name == 'nt':
            path_dest_log = join(getenv('APPDATA'), '.log', 'discord.log')
        else:
            path_dest_log = join(home, '.log', 'discord.log')
        bot.setup_logging(config['logging_level'], path_dest_log)
        bot.db = db_pool
        await bot.start(config['token'])

def run():
    asyncio.run(main())

if __name__ == '__main__':
    run()
