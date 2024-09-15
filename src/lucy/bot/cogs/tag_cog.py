''' tag_cog.py
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
import discord
from discord.ext import commands
from discord import app_commands
from typing import Annotated
import aiomysql

class TagCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.in_progress_tags = {}
        self.pool = bot.db_pool  # Ensure the bot has a db_pool attribute

    async def add_tag(self, ctx, name, content=None, attachment_url=None):
        query = """INSERT INTO tags (name, location_id, content, attachment_url, owner_id) 
                   VALUES (%s, %s, %s, %s, %s)"""
        async with self.pool.acquire() as conn:
            async with conn.cursor() as cursor:
                await cursor.execute(query, (name, ctx.guild.id, content, attachment_url, ctx.author.id))
                await conn.commit()

    async def get_tag(self, guild_id, name):
        query = """SELECT * FROM tags WHERE location_id=%s AND LOWER(name)=%s"""
        async with self.pool.acquire() as conn:
            async with conn.cursor(aiomysql.DictCursor) as cursor:
                await cursor.execute(query, (guild_id, name.lower()))
                tag = await cursor.fetchone()
                if tag:
                    return tag
                raise RuntimeError(f'Tag "{name}" not found.')

    async def update_tag(self, name, content=None, attachment_url=None, guild_id=None, owner_id=None):
        query = """UPDATE tags SET content = %s, attachment_url = %s 
                   WHERE name = %s AND location_id = %s AND owner_id = %s"""
        async with self.pool.acquire() as conn:
            async with conn.cursor() as cursor:
                result = await cursor.execute(query, (content, attachment_url, name, guild_id, owner_id))
                await conn.commit()
                return result

    async def delete_tag(self, name, guild_id, owner_id):
        query = """DELETE FROM tags WHERE name = %s AND location_id = %s AND owner_id = %s"""
        async with self.pool.acquire() as conn:
            async with conn.cursor() as cursor:
                result = await cursor.execute(query, (name, guild_id, owner_id))
                await conn.commit()
                return result

    @commands.hybrid_group(fallback='get')
    @commands.guild_only()
    @app_commands.guild_only()
    @app_commands.describe(name='The tag to retrieve')
    async def tag(self, ctx: commands.Context, *, name: Annotated[str, str]):
        """Allows you to retrieve a tag."""
        try:
            tag = await self.get_tag(ctx.guild.id, name)
        except RuntimeError as e:
            return await ctx.send(str(e))
    
        if tag['attachment_url']:
            await ctx.send(tag['attachment_url'])
        else:
            await ctx.send(tag['content'], reference=ctx.message.reference)
    
        query = "UPDATE tags SET uses = uses + 1 WHERE name = %s AND location_id=%s;"
        async with self.pool.acquire() as conn:
            async with conn.cursor() as cursor:
                await cursor.execute(query, (tag['name'], ctx.guild.id))

    @tag.command(aliases=['add'])
    @commands.guild_only()
    async def create(self, ctx: commands.Context, name: Annotated[str, str], *, content: Annotated[str, str] = None):
        """Creates a new tag with optional content or attachment."""
        if not ctx.message.attachments and content is None:
            return await ctx.send('You must provide either content or an attachment.')

        attachment_url = None
        if ctx.message.attachments:
            attachment_url = ctx.message.attachments[0].url

        await self.add_tag(ctx, name, content, attachment_url)
        await ctx.send(f'Tag "{name}" successfully created.')

    @tag.command()
    @commands.guild_only()
    async def edit(self, ctx: commands.Context, name: Annotated[str, str], *, content: Annotated[str, str] = None):
        """Edits an existing tag with optional content or attachment."""
        attachment_url = None
        if ctx.message.attachments:
            attachment_url = ctx.message.attachments[0].url

        result = await self.update_tag(name, content, attachment_url, ctx.guild.id, ctx.author.id)
        if result == 0:
            await ctx.send('Tag not found or you do not own this tag.')
        else:
            await ctx.send(f'Tag "{name}" successfully updated.')

    @tag.command()
    @commands.guild_only()
    async def remove(self, ctx: commands.Context, name: Annotated[str, str]):
        """Removes a tag owned by you."""
        result = await self.delete_tag(name, ctx.guild.id, ctx.author.id)
        if result == 0:
            await ctx.send('Tag not found or you do not own this tag.')
        else:
            await ctx.send(f'Tag "{name}" successfully removed.')

    @tag.command()
    @commands.guild_only()
    async def tags(self, ctx: commands.Context):
        """Lists all tags you own."""
        query = """SELECT name, content, attachment_url FROM tags WHERE location_id=%s AND owner_id=%s"""
        async with self.pool.acquire() as conn:
            async with conn.cursor() as cursor:
                await cursor.execute(query, (ctx.guild.id, ctx.author.id))
                tags = await cursor.fetchall()

        if not tags:
            await ctx.send('You have no tags.')
        else:
            tag_list = '\n'.join(f'{tag["name"]}: {tag["content"] if tag["content"] else tag["attachment_url"]}' for tag in tags)
            await ctx.send(f'Your tags:\n{tag_list}')

async def setup(bot):
    await bot.add_cog(TagCog(bot))
