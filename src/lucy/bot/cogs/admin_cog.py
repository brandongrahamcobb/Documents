''' admin_cog.py
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
from collections import defaultdict
from datetime import datetime, timedelta
from discord.ext import commands, tasks
from typing import Literal, Optional
from bot.main import Lucy

import asyncio
import bot.utils.helpers as lucy
import discord
import time
import json

CHANNEL_ID = 937463813626822656
GUILD_ID = 777341210871726090

def is_owner():
    async def predicate(ctx):
        return ctx.guild is not None and (ctx.guild.owner_id == ctx.author.id or ctx.author.id == 154749533429956608)
    return commands.check(predicate)

class AdminCog(commands.Cog):

    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.user_command_messages = {}
        self.persistent_users = lucy.path_users_yaml
        self.check_roles.start()  # Start the periodic task

    @commands.Cog.listener()
    async def on_member_join(self, member):
        now = datetime.now()
        user_data = {"joined_at": now.isoformat()}
        with open(self.persistent_users, 'a') as f:
            f.write(f"{member.id}:{json.dumps(user_data)}\n")

    @tasks.loop(minutes=5)
    async def check_roles(self):
        now = datetime.now()
        to_remove = []
        with open(self.persistent_users, 'r') as f:
            lines = f.readlines()
        for line in lines:
            user_id, data = line.split(':', 1)
            user_id = int(user_id)
            user_info = json.loads(data)
            joined_at = datetime.fromisoformat(user_info['joined_at'])
            time_difference = now - joined_at
            guild = await self.bot.fetch_guild(777341210871726090)
            member = guild.get_member(user_id)
            if member == None:
                break
            else:
                if member.roles > 1:
                    nerd_role = guild.get_role(1277654132177768609)
                    if nerd_role and nerd_role not in member.roles:
                        await member.add_roles(nerd_role)
                    to_remove.append(line)
                elif time_difference > timedelta(hours=24):
                    self.kick_function(member)
                    to_remove.append(line)
        with open(self.persistent_users, 'w') as f:
            for line in lines:
                if line not in to_remove:
                    f.write(line)

    def kick_function(self, member):
        print(f"{member.name} did not qualify for roles in 24 hours.")

    @check_roles.before_loop
    async def before_check_roles(self):
        await self.bot.wait_until_ready()

    @commands.Cog.listener()
    async def on_message_edit(self, before, after):
        if before.content != after.content:
            ctx = await self.bot.get_context(after)
            if ctx.command:
                await self.bot.invoke(ctx)

    @commands.Cog.listener()
    async def on_message(self, message: discord.Message) -> None:
        if message.author == self.bot.user:
            return
        # WIPE
        if message.author != self.bot.user and message.content.startswith('!'):
            if message.author.id not in self.user_command_messages:
                self.user_command_messages[message.author.id] = []
            self.user_command_messages[message.author.id].append(message.id)

    @commands.Cog.listener()
    async def on_ready(self):
        bot_user = self.bot.user
        bot_name = bot_user.name
        bot_id = bot_user.id
        guild_count = len(self.bot.guilds)
        info = (
            f'\n=============================\n'
            f'bot Name: {bot_name}\n'
            f'bot ID: {bot_id}\n'
            f'Connected Guilds: {guild_count}\n'
            f'============================='
        )
        guild_info = '\n'.join(
            [f'- {guild.name} (ID: {guild.id})' for guild in self.bot.guilds]
        )
        stats_message = f'{info}\n\nGuilds:\n{guild_info}'
        print(stats_message)
        self.delete_old_messages.start()

    @tasks.loop(hours=24)
    async def delete_old_messages(self):
        guild = await self.bot.fetch_guild(GUILD_ID)
        channel = guild.get_channel(CHANNEL_ID)
        if channel is not None:
            async for message in channel.history(limit=100):
                if (discord.utils.utcnow() - message.created_at).days >= 1:
                    await message.delete()
                    print(f'Deleted message from {message.author}: {message.content}')

    async def purge_messages(self, ctx, limit, check=None):
        deleted = 0
        async for message in ctx.channel.history(limit=limit):
            if check is None or check(message):
                await message.delete()
                deleted += 1
        return deleted

    @commands.command(name='load', hidden=True)
    @commands.is_owner()
    async def load(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.load_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.hybrid_command()
    async def reload(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.reload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.command(name='sync', hidden=True)
    @is_owner()
    async def sync(self, ctx: commands.Context, guilds: commands.Greedy[discord.Object], spec: Optional[Literal['~', '*', '^']] = None) -> None:
        if not guilds:
            if spec == '~':
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == '*':
                ctx.bot.tree.copy_global_to(guild=ctx.guild)
                synced = await ctx.bot.tree.sync(guild=ctx.guild)
            elif spec == '^':
                ctx.bot.tree.clear_commands(guild=ctx.guild)
                await ctx.bot.tree.sync(guild=ctx.guild)
                synced = []
            else:
                synced = await ctx.bot.tree.sync()
            await ctx.send(
                f'Synced {len(synced)} commands {'globally' if spec is None else 'to the current guild.'}'
            )
            return
        ret = 0
        for guild in guilds:
            try:
                await ctx.bot.tree.sync(guild=guild)
            except discord.HTTPException:
                pass
            else:
                ret += 1
        await ctx.send(f'Synced the tree to {ret}/{len(guilds)}.')

    @commands.command(name='unload', hidden=True)
    @commands.is_owner()
    async def unload(self, ctx: commands.Context, *, module: str):
        try:
            await ctx.bot.unload_extension(module)
        except commands.ExtensionError as e:
            await ctx.send(f'{e.__class__.__name__}: {e}')
        else:
            await ctx.send('\N{OK HAND SIGN}')

    @commands.command(name='wipe', hidden=True)
    @is_owner()
    async def wipe(self, ctx, option=None, limit: int = 100):
        if limit <= 0 or limit > 100:
            await ctx.send('Please provide a limit between 1 and 100.')
            return
        if option == 'bot':
            def is_bot_message(message):
                return message.author == self.bot.user
            deleted = await self.purge_messages(ctx, limit, is_bot_message)
            await ctx.send(f'Deleted {deleted} bot messages.')
        elif option == 'all':
            deleted = await ctx.channel.purge(limit=limit)
            await ctx.send(f'Deleted all messages.')
        elif option == 'user':
            def is_user_message(message):
                return message.author == ctx.author
            deleted = await self.purge_messages(ctx, limit, is_user_message)
            await ctx.send(f'Deleted {deleted} of your messages.')
        elif option == 'commands':
            if ctx.author.id in self.user_command_messages:
                message_ids = self.user_command_messages[ctx.author.id]
                deleted = 0
                for message_id in message_ids:
                    try:
                        message = await ctx.channel.fetch_message(message_id)
                        await message.delete()
                        deleted += 1
                    except discord.NotFound:
                        continue
                del self.user_command_messages[ctx.author.id]
                await ctx.send(f'Deleted {deleted} of your commands.')
            else:
                await ctx.send('No commands found to delete.')
        else:
            await ctx.send('Invalid option. Use `bot`, `all`, `user`, or `commands`.')

async def setup(bot: commands.bot):
    await bot.add_cog(AdminCog(bot))
